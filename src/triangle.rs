use std::collections::HashSet;

use super::vertex::Vertex;
use super::edge::Edge;

#[derive(Hash, Eq, PartialEq, Debug)]
struct Triangle {
	edge: &Edge,
}

impl Index for Triangle {
	fn index(&mut self, index: u32) -> &Vertex {
		match index {
			0 => self.edge.start(),
			1 => self.edge.end(),
			2 => self.edge.next().end(),
			_ => panic!("Unknown point {} (triangles only have three points dummy)", index),
		}
	}
}

impl Triangle {
	pub fn new(edge: &Edge) -> Triangle { Triangle { edge: edge } }
	/*
	pub fn neighboors(&self) -> Vec<&Triangles> {
		let mut neighboors: Vec<&Triangles> = Vec::new();
		
		let mut current_edge: Option<&Edge> = self.edge.opposite();
		// start searching in one direction
		while current_edge != self.edge  {
			current_edge.next().opposite();
			
		}
		// search in the other direction if the first direction was null
		if current_edge == None {
			current_edge = self.edge.next().opposite();
			while current_edge != self.edge  {
				current_edge = current_edge.next().opposite();
				
			}
		}
		
		neighboors;
	}*/
}
// edge()