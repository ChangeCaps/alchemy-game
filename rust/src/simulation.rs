use std::{
    hash::{Hash, Hasher},
    mem,
};

use godot::prelude::*;
use rayon::prelude::*;
use real_consts::PI;

macro_rules! time {
    ($name:expr, $block:expr) => {{
        let start = std::time::Instant::now();
        let result = $block;
        let elapsed = start.elapsed();
        godot_print!("{}: {:?}", $name, elapsed);
        result
    }};
}

#[derive(GodotClass, Default)]
#[class(init)]
pub struct Simulation {
    #[var]
    #[init(default = 1.0)]
    pub pressure_multiplier: f32,

    #[var]
    #[init(default = 1.0)]
    pub target_density: f32,

    #[var]
    #[init(default = 1.0)]
    pub smoothing_radius: f32,

    #[var]
    #[init(default = 5.0)]
    pub boundary_radius: f32,

    #[var]
    #[init(default = 0.0)]
    pub boundary_density: f32,

    atoms: Vec<AtomData>,
    particles: Vec<Particle>,
    densities: Vec<f32>,
    pressures: Vec<Vector2>,

    spatial_lut: Vec<u16>,
    lut_indices: Vec<u16>,
}

#[godot_api]
impl Simulation {
    #[func]
    fn set_atoms(&mut self, atoms: Array<Gd<Atom>>) {
        self.atoms = atoms.iter_shared().map(|a| a.bind().data()).collect();
    }

    #[func]
    fn get_atoms(&self) -> Array<Gd<Atom>> {
        self.atoms
            .iter()
            .map(|a| {
                Gd::from_init_fn(|base| Atom {
                    bonds: a.bonds,
                    mass: a.mass,
                    charge: a.charge,
                    base,
                })
            })
            .collect()
    }

    #[func]
    fn add_particle(&mut self, particle: Gd<Particle>) {
        self.particles.push(*particle.bind());
    }

    #[func]
    fn particle_count(&self) -> i64 {
        self.particles.len() as i64
    }

    #[func]
    fn particle_position(&self, index: i64) -> Vector2 {
        self.particles[index as usize].position
    }

    #[func]
    fn average_density(&self) -> f32 {
        self.densities.iter().sum::<f32>() / self.densities.len() as f32
    }

    #[func]
    fn step(&mut self, dt: f32) {
        time!("update_spatial_lut", self.update_spatial_lut(dt));
        time!("update_densities", self.update_densities(dt));
        time!("update_pressures", self.update_pressures(dt));
        time!("update_velocities", self.update_velocities(dt));
        time!("update_positions", self.update_positions(dt));
        time!("resolve_boundary", self.resolve_boundary());
    }

    fn update_spatial_lut(&mut self, dt: f32) {
        self.spatial_lut.resize(self.particles.len(), 0);
        let mut spatial_lut = mem::take(&mut self.spatial_lut);

        spatial_lut.iter_mut().enumerate().for_each(|(i, hash)| {
            let p = &self.particles[i];
            let point = p.position + p.velocity * dt;
            let cell = self.compute_cell(point);
            *hash = self.compute_cell_hash(cell);
        });

        self.spatial_lut = spatial_lut;
        self.spatial_lut.sort_unstable();

        self.lut_indices.resize(self.particles.len(), 0);

        if let Some(first) = self.spatial_lut.first() {
            self.lut_indices[*first as usize] = 0;
        }

        for i in 1..self.spatial_lut.len() {
            let hash = self.spatial_lut[i];
            let prev = self.spatial_lut[i - 1];

            if hash != prev {
                self.lut_indices[hash as usize] = i as u16;
            }
        }

        godot_print!("spatial {:?}", self.spatial_lut);
        godot_print!("indices {:?}", self.lut_indices);
    }

    fn update_densities(&mut self, dt: f32) {
        self.densities.resize(self.particles.len(), 0.0);
        let mut densities = mem::take(&mut self.densities);

        self.particles
            .iter()
            .zip(densities.iter_mut())
            .for_each(|(p, d)| {
                let point = p.position + p.velocity * dt;

                *d = self.compute_density(point, dt);
            });

        self.densities = densities;
    }

    fn update_pressures(&mut self, dt: f32) {
        self.pressures.resize(self.particles.len(), Vector2::ZERO);
        let mut pressures = mem::take(&mut self.pressures);

        pressures.iter_mut().enumerate().for_each(|(i, pr)| {
            *pr = self.compute_pressure(i, dt);
        });

        self.pressures = pressures;
    }

    fn update_velocities(&mut self, dt: f32) {
        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, p)| {
                let acceleration = self.pressures[i] / self.densities[i];

                p.velocity += acceleration * dt;
            });
    }

    fn update_positions(&mut self, dt: f32) {
        self.particles.par_iter_mut().for_each(|p| {
            p.position += p.velocity * dt;
        });
    }

    fn resolve_boundary(&mut self) {
        self.particles.par_iter_mut().for_each(|p| {
            if p.position.length() > self.boundary_radius {
                p.position = p.position.normalized() * self.boundary_radius;
                p.velocity = -p.velocity.reflect(-Vector2::normalized(p.position));
                p.velocity *= 0.5;
            }
        });
    }

    fn compute_pressure(&self, index: usize, dt: f32) -> Vector2 {
        let a = self.particles[index];
        let a_point = a.position + a.velocity * dt;

        self.particles_in_radius(a_point)
            .map(|(i, b)| {
                let b_point = b.position + b.velocity * dt;

                let distance = a_point.distance_to(b_point);

                if distance < 1e-6 {
                    return Vector2::ZERO;
                }

                let atom = &self.atoms[b.atom as usize];
                let direction = Vector2::normalized(a_point - b_point);

                let gradient = Self::smooth_kernel_gradient(self.smoothing_radius, distance);

                let a_density = self.densities[index];
                let b_density = self.densities[i];

                let a_pressure = self.density_to_pressure(a_density);
                let b_pressure = self.density_to_pressure(b_density);

                let pressure = (a_pressure + b_pressure) / 2.0;

                -direction * pressure * gradient * atom.mass / b_density
            })
            .sum()
    }

    fn compute_density(&self, point: Vector2, dt: f32) -> f32 {
        let boundary_distance = f32::max(0.0, self.smoothing_radius - point.length());
        let boundary_density = Self::smooth_kernel(self.smoothing_radius, boundary_distance);
        let boundary_density = boundary_density * self.boundary_density;

        self.particles_in_radius(point)
            .map(|(_, p)| {
                let p_point = p.position + p.velocity * dt;

                let atom = &self.atoms[p.atom as usize];
                let distance = p_point.distance_to(point);

                Self::smooth_kernel(self.smoothing_radius, distance) * atom.mass
            })
            .sum::<f32>()
            + boundary_density
    }

    fn particles_in_radius(&self, point: Vector2) -> impl Iterator<Item = (usize, &Particle)> {
        const NEIGHBORS: [Vector2i; 9] = [
            Vector2i::new(-1, -1),
            Vector2i::new(-1, 0),
            Vector2i::new(-1, 1),
            Vector2i::new(0, -1),
            Vector2i::new(0, 0),
            Vector2i::new(0, 1),
            Vector2i::new(1, -1),
            Vector2i::new(1, 0),
            Vector2i::new(1, 1),
        ];

        let cell = self.compute_cell(point);

        NEIGHBORS.into_iter().flat_map(move |offset| {
            let cell = cell + offset;
            let hash = self.compute_cell_hash(cell);
            let index = self.lut_indices[hash as usize] as usize;

            self.spatial_lut[index..]
                .iter()
                .enumerate()
                .take_while(move |(_, lut_hash)| **lut_hash == hash)
                .map(move |(i, _)| {
                    let particle = &self.particles[index + i];
                    (index + i, particle)
                })
        })
    }

    fn density_to_pressure(&self, density: f32) -> f32 {
        let delta_density = density - self.target_density;
        self.pressure_multiplier * delta_density
    }

    fn compute_cell(&self, point: Vector2) -> Vector2i {
        let cell = point / self.smoothing_radius;

        Vector2i::new(cell.x.floor() as i32, cell.y.floor() as i32)
    }

    fn compute_cell_hash(&self, cell: Vector2i) -> u16 {
        let mut hasher = seahash::SeaHasher::new();

        cell.hash(&mut hasher);

        (hasher.finish() as usize % self.particles.len()) as u16
    }

    fn smooth_kernel(r: f32, d: f32) -> f32 {
        if d >= r {
            return 0.0;
        }

        let volume = (PI * r.powi(4)) / 6.0;
        (r - d) * (r - d) / volume
    }

    fn smooth_kernel_gradient(r: f32, d: f32) -> f32 {
        const EPSILON: f32 = 1e-6;

        let delta = Self::smooth_kernel(r, d + EPSILON) - Self::smooth_kernel(r, d - EPSILON);
        delta / (2.0 * EPSILON)
    }
}

#[derive(Clone, Copy, Debug, Default, GodotClass)]
#[class(init)]
pub struct Particle {
    #[export]
    pub position: Vector2,

    #[export]
    pub velocity: Vector2,

    pub bonds: [u16; 4],

    #[export]
    pub atom: u16,
}

impl Particle {
    pub fn bond_count(&self) -> u8 {
        let [a, b, c, d] = self.bonds;

        let a = (a == u16::MAX) as u8;
        let b = (b == u16::MAX) as u8;
        let c = (c == u16::MAX) as u8;
        let d = (d == u16::MAX) as u8;

        a + b + c + d
    }

    pub fn bonds(&self) -> impl Iterator<Item = u16> {
        self.bonds.into_iter().filter(|&b| b != u16::MAX)
    }
}

struct AtomData {
    bonds: u8,
    mass: f32,
    charge: f32,
}

#[derive(GodotClass)]
#[class(base = Resource)]
pub struct Atom {
    #[export]
    pub bonds: u8,

    #[export]
    pub mass: f32,

    #[export]
    pub charge: f32,

    base: Base<Resource>,
}

#[godot_api]
impl IResource for Atom {
    fn init(base: Base<Resource>) -> Self {
        Self {
            bonds: 0,
            mass: 1.0,
            charge: 0.0,
            base,
        }
    }
}

impl Atom {
    fn data(&self) -> AtomData {
        AtomData {
            bonds: self.bonds,
            mass: self.mass,
            charge: self.charge,
        }
    }
}
