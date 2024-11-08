use godot::prelude::*;

struct RustExtension;

#[gdextension]
unsafe impl ExtensionLibrary for RustExtension {}

#[derive(GodotClass)]
#[class(base = Node3D)]
struct Test {
    #[export]
    speed: f64,

    base: Base<Node3D>,
}

#[godot_api]
impl INode3D for Test {
    fn init(base: Base<Node3D>) -> Self {
        Test { speed: 5.0, base }
    }
}
