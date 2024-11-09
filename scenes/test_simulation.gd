extends Node3D

@export var pressure_multiplier: float = 1.0
@export var target_density: float = 1.0
@export var smoothing_radius: float = 1.0
@export var boundary_radius: float = 5.0
@export var boundary_density: float = 0.0
@export var atoms: Array[Atom]

var simulation: Simulation = Simulation.new()

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
    simulation.set_atoms(atoms)

    for i in range(10):
        var particle = Particle.new()
        particle.position = Vector2(randf(), randf())

        simulation.add_particle(particle)


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
    simulation.pressure_multiplier = pressure_multiplier
    simulation.target_density = target_density
    simulation.smoothing_radius = smoothing_radius
    simulation.boundary_radius = boundary_radius
    simulation.boundary_density = boundary_density

    $MultiMeshInstance3D.multimesh.instance_count = simulation.particle_count()

    simulation.step(delta)

    for i in range(simulation.particle_count()):
        var p = simulation.particle_position(i)
        var t = Transform3D(Basis(), Vector3(p.x, 0.0, p.y))

        $MultiMeshInstance3D.multimesh.set_instance_transform(i, t)
