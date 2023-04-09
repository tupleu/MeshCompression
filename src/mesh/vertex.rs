use std::fmt;

#[derive(Hash, Eq, PartialEq, Debug, Clone, Default)]
pub struct Vertex {
	x: i32,
	y: i32,
}

impl fmt::Display for Vertex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Vertex ({}, {})", self.x(), self.y())
    }
}

impl Vertex {
	pub fn new(x: i32, y: i32) -> Vertex { Vertex { x: x, y: y } }
	
	pub fn x(&self) -> i32 { self.x }
	pub fn y(&self) -> i32 { self.y }
	
	pub fn pos(&self) -> (i32, i32) { (self.x, self.y) }
}