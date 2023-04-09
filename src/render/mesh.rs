use std::collections::HashMap;

pub mod vertex;
mod edge;
mod triangle;
use self::vertex::{Point, Vertex};
use self::edge::Edge;
//use self::triangle::Triangle;

pub struct Mesh {
	verticies: Vec<Vertex>,
	edges: Vec<Edge>,
	triangles: Vec<u16>,
	vertex_map: HashMap<Point, u16>,
	edge_map: HashMap<(u16, u16), u16>,
}

impl Mesh {
	pub fn new() -> Mesh { Mesh { verticies: vec![], edges: vec![], triangles: vec![], vertex_map: HashMap::new(), edge_map: HashMap::new() } }
	
	pub fn from_image(width: usize, height: usize) -> Mesh {
		if width < 2 || height < 2 {
			return Mesh::new();
		}
		
		let mut mesh = Mesh::new();
		
		let mut vs = vec![vec![0; height+1]; width + 1];
		// make the verticies
		for i in 0..width {
			for j in 0..height {
				vs[i][j] = mesh.verticies.len();
				mesh.get_or_add_vertex(i as f32, j as f32, Some([1.0, 0.0, 0.0]));
			}
		}
		
		// create the edges and faces
		for i in 0..width {			
			let i = i as f32;
			for j in 0..height {
				let j = j as f32;
				mesh.add_triangle((i,j), (i,j+1.0), (i+1.0,j));
				mesh.add_triangle((i,j+1.0), (i+1.0,j), (i+1.0,j+1.0));
			}
		}
		//mesh.add_triangle([vs[0][0], vs[1][0], vs[0][1]]);
		
		mesh
	}
	
	pub fn verticies(&self) -> &Vec<Vertex> {
		return &self.verticies;
	}
	
	pub fn indicies(&self) -> &Vec<u16> {
		return &self.triangles;
	}
	
	pub fn get_or_add_vertex(&mut self, x: f32, y: f32, color: Option<[f32; 3]>) -> u16 {
		if let Some(index) = self.vertex_map.get(&Point::new(x, y)) {
			return *index;
		}
		let index = self.verticies.len() as u16;
		self.verticies.push( Vertex::new(x, y, color.unwrap_or([0.0, 0.0, 0.0])) );
		self.vertex_map.insert(Point::new(x, y), index);
		index
	}
	
	pub fn get_vertex(&self, index: u16) -> &Vertex {
		&self.verticies[index as usize]
	}
	
	pub fn get_edge(&self, index: u16) -> &Edge {
		&self.edges[index as usize]
	}
	
	pub fn get_edge_mut(&mut self, index: u16) -> &mut Edge {
		&mut self.edges[index as usize]
	}
	
	pub fn has_vertex(&self, v: (f32, f32)) -> i32 {
		match self.vertex_map.get(&Point::new(v.0, v.1)) {
			Some(index) => *index as i32,
			None => -1,
		}
	}
	
	pub fn has_edge(&self, v_start: (f32, f32), v_end: (f32, f32)) -> i32 {
		let v1_index = self.has_vertex(v_start);
		let v2_index = self.has_vertex(v_end);
		if v1_index >= 0 && v2_index >= 0 {
			return match self.edge_map.get(&(v1_index as u16, v2_index as u16)) {
				Some(index) => *index as i32,
				None => match self.edge_map.get(&(v2_index as u16, v1_index as u16)) {
							Some(index) => *index as i32,
							None => -1,
						},
			};
		}
		-1
	}
	
	fn add_triangle(&mut self, mut v1: (f32, f32), mut v2: (f32, f32), mut v3: (f32, f32) ) {
		let ei = [self.has_edge(v1, v2), self.has_edge(v2, v3), self.has_edge(v3, v1)];
		
		if ei[0] >= 0 {
			let e_temp = self.get_edge(ei[0] as u16);
			v1 = self.get_vertex(e_temp.get_end()).pos();
			v2 = self.get_vertex(e_temp.get_start()).pos();
		} else if ei[1] >= 0 {
			let e_temp = self.get_edge(ei[1] as u16);
			v2 = self.get_vertex(e_temp.get_end()).pos();
			v3 = self.get_vertex(e_temp.get_start()).pos();
		} else if ei[2] >= 0 {
			let e_temp = self.get_edge(ei[2] as u16);
			v3 = self.get_vertex(e_temp.get_end()).pos();
			v1 = self.get_vertex(e_temp.get_start()).pos();
		}
		
		let vertecies = [self.get_or_add_vertex(v1.0, v1.1, None),
						 self.get_or_add_vertex(v2.0, v2.1, None),
						 self.get_or_add_vertex(v3.0, v3.1, None)];
						 
		let edges_len = self.edges.len();
		
		for index in 0..3 {
			let index_clamped = (index + 1) % 3;
			let next_edge_index = index_clamped + edges_len;
			let v_start = vertecies[index] as u16;
			let v_end = vertecies[index_clamped] as u16;
			let mut edge = Edge::new(v_start, v_end, next_edge_index as u16);
			edge.set_opposite(ei[index]);
			if ei[index] >= 0 {
				self.get_edge_mut(ei[index] as u16).set_opposite((edges_len + index) as i32);
			}
			
			
			self.edges.push(edge);
			self.edge_map.insert((v_start, v_end), (edges_len + index) as u16);
		}
		self.triangles.push(vertecies[0]);
		self.triangles.push(vertecies[1]);
		self.triangles.push(vertecies[2]);
	}
	/*
	fn add_triangle_to_edge(&mut self, edge_index: usize, vertex: (i32, i32)) {
		let vi1 = self.get_edge(edge_index).get_end();
		let vi2 = self.get_edge(edge_index).get_start();
		
		self.add_triangle((vi1.x(), vi1.y()), (vi2.x(), vi2.y()), vertex);
	}
	*/
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