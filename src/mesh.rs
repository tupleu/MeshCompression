use std::collections::HashMap;
use std::collections::LinkedList;

mod vertex;
mod edge;
mod triangle;
use self::vertex::Vertex;
use self::edge::Edge;
use self::triangle::Triangle;

pub struct Mesh {
	verticies: Vec<Vertex>,
	edges: Vec<Edge>,
	triangles: Vec<Triangle>,
	vertex_map: HashMap<(i32, i32), usize>,
	edge_map: HashMap<(usize, usize), usize>,
}

impl Mesh {
	pub fn new() -> Mesh { Mesh { verticies: vec![], edges: vec![], triangles: vec![], vertex_map: HashMap::new(), edge_map: HashMap::new() } }
	
	pub fn from_image(width: u32, height: u32) -> Mesh { 
		Mesh { verticies: vec![], edges: vec![], triangles: vec![], vertex_map: HashMap::new(), edge_map: HashMap::new() }
	}
	
	pub fn get_or_add_vertex(&mut self, x: i32, y: i32) -> &Vertex {
		if let Some(index) = self.vertex_map.get(&(x, y)) {
			return &self.get_vertex(*index);
		}
		let index = self.verticies.len();
		self.verticies.push( Vertex::new(x, y) );
		self.vertex_map.insert((x, y), index);
		&self.verticies[index]
	}
	
	pub fn get_vertex(&self, index: usize) -> &Vertex {
		&self.verticies[index]
	}
	
	pub fn get_edge(&self, index: usize) -> &Edge {
		&self.edges[index]
	}
	
	fn add_triangle(&mut self, vertecies: [usize; 3]) {
		let edges_len = self.edges.len();
		for index in 0..3 {
			let index_clamped = (index + 1) % 3;
			let next_edge_index = index_clamped + edges_len;
			let v_start = vertecies[index];
			let v_end = vertecies[index_clamped];
			let edge = Edge::new(v_start, v_end, next_edge_index);
			self.edges.push(edge);
			self.edge_map.insert((v_start, v_end), edges_len + index);
		}
		let tri = Triangle::new(edges_len);
		self.triangles.push(tri);
	}
	
	fn add_triangle_to_edge(&mut self, edge_index: usize, vertex_index: usize) {
		let vi1 = self.get_edge(edge_index).get_end();
		let vi2 = self.get_edge(edge_index).get_start();
		
		self.add_triangle([vi1, vi2, vertex_index]);
	}
	
	//pub fn length(&self, edge: usize) -> f64 { ( ( i32::pow(self.get() - self.start.x(), 2) + i32::pow(self.end.y() - self.start.y(), 2) ) as f64 ).sqrt() }
	/*
	pub fn init_from_image(&mut self, width: u32, height: u32) { 
		let mut prev_e: LinkedList<Edge> = LinkedList::new();
		let mut prev_v: LinkedList<Vertex> = LinkedList::new();
		
		
		for x in 0..width {
			for y in 0..height {
				if prev_v.len() == 2 {
					let v1 = prev_v.back().expect("rust sucks 1");
					let v2 = prev_v.front().expect("rust sucks 2");
					let e = Edge::new(v1, v2);
					prev_e.push_front(e);
					//prev_e.push_front(Rc::new(Edge::new(prev_v.back().expect("Vertex not in set"), prev_v.front().expect("Vertex not in set"))));
					
					//v.insert(prev_v.pop_back().expect("rust sucks"));
					self.verticies.insert(prev_v.pop_back().expect("rust sucks"));
				}
				let v = Vertex::new(x as i32, y as i32);
				prev_v.push_front(v);
			}
		}
		
		
		
	}
	
	*/
}
