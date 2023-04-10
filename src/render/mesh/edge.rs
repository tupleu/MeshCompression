use super::vertex::Vertex;

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
pub struct Edge {
	start: u32,
	end: u32,
	next: u32, 
	opposite: Option<u32>,
}

impl Edge {
	
	pub fn new(start: u32, end: u32, next: u32) -> Edge { Edge { start: start, end: end, next: next, opposite: None } }
	
	pub fn get_start(&self) -> u32 { self.start }
	pub fn get_end(&self) -> u32 { self.end }
	
	pub fn set_opposite(&mut self, edge: i32) {
		if edge >= 0 {
			self.opposite = Some(edge as u32);
		}
	}
	
	//pub fn start(&self) -> &Vertex { self.start }
	//pub fn end(&self) -> &Vertex { self.end }
	//pub fn triangle(&self) -> &Triangle { self.triangle }
	
	//pub fn next(&self) -> Option<&Edge> { self.next }
	//pub fn prev(&self) -> Option<&Edge> { self.next?.next() }
	//pub fn opposite(&self) -> Option<&Edge> { self.opposite }
	
	//pub fn length(&self) -> f64 { ( ( i32::pow(self.end.x() - self.start.x(), 2) + i32::pow(self.end.y() - self.start.y(), 2) ) as f64 ).sqrt() }
	
}

// start(), end(), next(), opposite(), lenght(), triangle()
// Nuclear Naval Labratories

// Global Foundry, team 











