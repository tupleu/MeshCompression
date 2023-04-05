mod vertex;
mod edge;

use vertex::Vertex;
use edge::Edge;

fn main() {
    println!("Hello, world!");
	let v = Vertex::new(1, 1);
	println!("v: ({}, {})", v.x(), v.y());
	println!("{}", v);
}
