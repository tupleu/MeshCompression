
static EPSILON: f32 = f32::EPSILON * 4.0;
static MULT: f64 = 0x5d9f3b as f64;
static MULT1: i64 = 0x9de1f3;

fn float_hash(num: f32) -> i64 {
	(num as f64 * MULT) as i64 + MULT1
}

fn float_tuple_hash(num: (f32, f32)) -> i64 {
	(float_hash(num.0) * MULT1) + float_hash(num.1)
}

fn float_tuple_eq(a: (f32, f32), b: (f32, f32)) -> bool {
	(a.0 - b.0).abs() < EPSILON && (a.1 - b.1).abs() < EPSILON
}

#[derive(Copy, Clone, Debug, Default)]
pub struct Point {
	pub x: f32,
	pub y: f32,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Copy)]
pub struct VertexIndex { pub index: u32 }
impl VertexIndex { pub fn new(index: u32) -> VertexIndex { VertexIndex { index: index } } }

impl PartialEq for Point {
    fn eq(&self, other: &Point) -> bool {
        if 		(self.x - other.x).abs() < EPSILON && 
				(self.y - other.y).abs() < EPSILON		 {
			return true;
		}
		false
    }
}

impl Eq for Point {}

impl std::hash::Hash for Point {
    fn hash<H>(&self, state: &mut H)
    where
        H: std::hash::Hasher,
    {
		//let hash = (ival_x * MULT1  * MULT1) + (ival_y * MULT1) + zIntValue;
        state.write_i64(float_tuple_hash((self.x, self.y)));
        state.finish();
    }
}

impl Point {
	pub fn new(x: f32, y: f32) -> Point { Point { x: x, y: y } }
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub(crate) position: [f32; 3],
    pub(crate) color: [f32; 3],
}

impl PartialEq for Vertex {
    fn eq(&self, other: &Vertex) -> bool {
        if 		(self.position[0] - other.position[0]).abs() < EPSILON && 
				(self.position[1] - other.position[1]).abs() < EPSILON		 {
			return true;
		}
		false
    }
}

impl Eq for Vertex {}

impl std::hash::Hash for Vertex {
    fn hash<H>(&self, state: &mut H)
    where
        H: std::hash::Hasher,
    {
		//let hash = (ival_x * MULT1  * MULT1) + (ival_y * MULT1) + zIntValue;
        state.write_i64(float_tuple_hash((self.x(), self.y())));
        state.finish();
    }
}

impl std::fmt::Display for Vertex {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Vertex ({}, {})", self.x(), self.y())
    }
}

impl Vertex {
	pub fn new(x: f32, y: f32, color: [f32; 3]) -> Vertex { Vertex { position: [x, y, 0.0], color: color } }
	
	
	pub fn x(&self) -> f32 { self.position[0] }
	pub fn y(&self) -> f32 { self.position[1] }
	
	pub fn pos(&self) -> Point { Point::new(self.x(), self.y()) }
	pub fn color(&self) -> &[f32; 3] { &self.color }
	
	pub(crate) fn set_color(&mut self, color: [f32; 3]) { self.color = color }
}