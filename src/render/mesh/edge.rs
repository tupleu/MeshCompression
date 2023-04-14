use super::vertex::VertexIndex;

#[derive(Hash, Eq, PartialEq, Debug, Clone, Copy)]
pub struct EdgeIndex { pub index: u32 }
impl EdgeIndex { pub fn new(index: u32) -> EdgeIndex { EdgeIndex { index: index } } }


#[derive(Hash, Eq, PartialEq, Debug, Clone)]
pub struct Edge {
	start: VertexIndex,
	end: VertexIndex,
	next: EdgeIndex, 
	opposite: Option<EdgeIndex>,
}


impl Edge {
	
	pub fn new(start: VertexIndex, end: VertexIndex, next: EdgeIndex) -> Edge { Edge { start: start, end: end, next: next, opposite: None } }
	
	pub fn start(&self) -> VertexIndex { self.start }
	pub fn end(&self) -> VertexIndex { self.end }
	
	pub(crate) fn set_opposite(&mut self, edge: Option<EdgeIndex>) {
		self.opposite = edge;
	}
	
	pub fn opposite(&self) -> Option<EdgeIndex> {
		self.opposite
	}
	
	pub fn next(&self) -> EdgeIndex {
		self.next
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











