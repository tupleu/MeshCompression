mod vertex;

struct Mesh {
	verticies: HashSet<Vertex>,
	edges: HashSet<Edge>,
	triangles: HashSet<Triangle>,
}

impl Index for Mesh {
	fn index(&mut self, index: u32) -> &Vertex {
		match index {
			0 => 
			1 => 
			2 => 
			_ => panic!("Unknown point {} (triangles only have three points dummy)", index),
		}
	}
}

impl Mesh {
	pub fn new() -> Mesh { Vertex { verticies: HashSet<Vertex>::new(), edges: HashSet<Edge>::new(), triangles: HashSet<Triangle>::new() } }
	
	
	pub fn addTriangle(&self) -> i32 { self.x }
	pub fn remove(&self) -> i32 { self.y }
	
}
