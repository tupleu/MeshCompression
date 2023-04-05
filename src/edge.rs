use std::collections::HashSet;

use super::vertex::Vertex;
use super::triangle::Triangle;

#[derive(Hash, Eq, PartialEq, Debug)]
struct Edge {
	start: &Vertex,
	end: &Vertex,
	triangle: &Triangle,
	next: &Edge, 
	opposite: &Edge,
}

impl Edge {
	pub fn new(start: &Vertex, end: &Vertex) -> Edge { Edge { start: start, end: end } }
	
	pub fn start(&self) -> i32 { self.start }
	pub fn end(&self) -> i32 { self.end }
	pub fn triangle(&self) -> &Triangle { self.triangle }
	
	pub fn next(&self) -> &Edge { self.next }
	pub fn prev(&self) -> &Edge { self.next.next() }
	pub fn opposite(&self) -> Option<&Edge> { self.opposite }
	
	pub fn length(&self) -> f64 { ( ( i32::pow(self.end.x() - self.start.x(), 2) + i32::pow(self.end.y() - self.start.y(), 2) ) as f64 ).sqrt() }
}

// start(), end(), next(), opposite(), lenght(), triangle()
// Nuclear Naval Labratories

// Global Foundry, team 




































sample size
frequency
	- bl_service.py (line 58)
	- https://docs.python.org/3.5/library/struct.html#module-struct
	- bt-central.py (line 100, 47, 115)
	- settings.html (line 139ish)

Save 1 set of data for each time we sample
	- bt-central.py (line 233)
	- in a new function, we will grab the relevant data
		- throw away 0s
		- calculate avg, std deviation, average flourescence factor, (slope)
			- for the slope, we want a function that fits the expected V vs t graph
			- https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
		- also record the number of samples

study mode

ability to set a sample

throw 0s 

avg, std deviation, # of samples, ff, slope, 

std dev w/ 0s?


'''study mode''' (store in different database)

turn on and off saving all samples







