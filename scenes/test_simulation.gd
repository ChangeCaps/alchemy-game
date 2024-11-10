extends Node3D

@export var pressure_multiplier: float = 1.0
@export var target_density: float = 1.0
@export var smoothing_radius: float = 1.0
@export var boundary_radius: float = 5.0
@export var boundary_density: float = 0.0
@export var atoms: Array[Atom]

var simulation: Simulation = Simulation.new()

@onready var mesh: MultiMeshInstance3D = $MultiMeshInstance3D

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
    simulation.set_atoms(atoms)

    for i in range(5000):
        var particle = Particle.new()
        particle.position = Vector2(randf(), randf()) * 3.0 - Vector2(1.5, 1.5)

        simulation.add_particle(particle)


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
    simulation.pressure_multiplier = pressure_multiplier
    simulation.target_density = target_density
    simulation.smoothing_radius = smoothing_radius
    simulation.boundary_radius = boundary_radius

    mesh.multimesh.instance_count = simulation.particle_count()

    simulation.step(delta)

    for i in range(simulation.particle_count()):
        var p = simulation.particle_position(i)
        var t = Transform3D(Basis(), Vector3(p.x, 0.0, p.y))

        var density = simulation.particle_density(i)

        var color

        if density < target_density:
            color = Color(0.0, 0.0, target_density - density) / 5.0
        else:
            color = Color(density - target_density, 0.0, 0.0) / 5.0

        mesh.multimesh.set_instance_color(i, color) 
        mesh.multimesh.set_instance_transform(i, t)
