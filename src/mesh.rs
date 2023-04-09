use std::collections::HashMap;

mod vertex;
mod edge;
mod triangle;
use self::vertex::Vertex;
use self::edge::Edge;
use self::triangle::Triangle;

pub struct Mesh {
	verticies: Vec<Vertex>,
	edges: Vec<Edge>,
	pub triangles: Vec<Triangle>,
	vertex_map: HashMap<(i32, i32), usize>,
	edge_map: HashMap<(usize, usize), usize>,
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
				mesh.get_or_add_vertex(i as i32, j as i32);
			}
		}
		
		// create the edges and faces
		for i in 0..width {			
			let i = i as i32;
			for j in 0..height {
				let j = j as i32;
				mesh.add_triangle((i,j), (i,j+1), (i+1,j));
				mesh.add_triangle((i,j+1), (i+1,j), (i+1,j+1));
			}
		}
		//mesh.add_triangle([vs[0][0], vs[1][0], vs[0][1]]);
		
		mesh
	}
	
	pub fn get_or_add_vertex(&mut self, x: i32, y: i32) -> usize {
		if let Some(index) = self.vertex_map.get(&(x, y)) {
			return *index;
		}
		let index = self.verticies.len();
		self.verticies.push( Vertex::new(x, y) );
		self.vertex_map.insert((x, y), index);
		index
	}
	
	pub fn get_vertex(&self, index: usize) -> &Vertex {
		&self.verticies[index]
	}
	
	pub fn get_edge(&self, index: usize) -> &Edge {
		&self.edges[index]
	}
	
	pub fn get_edge_mut(&mut self, index: usize) -> &mut Edge {
		&mut self.edges[index]
	}
	
	pub fn has_vertex(&self, v: (i32, i32)) -> i32 {
		match self.vertex_map.get(&v) {
			Some(index) => *index as i32,
			None => -1,
		}
	}
	
	pub fn has_edge(&self, v_start: (i32, i32), v_end: (i32, i32)) -> i32 {
		let v1_index = self.has_vertex(v_start);
		let v2_index = self.has_vertex(v_end);
		if v1_index >= 0 && v2_index >= 0 {
			return match self.edge_map.get(&(v1_index as usize, v2_index as usize)) {
				Some(index) => *index as i32,
				None => match self.edge_map.get(&(v2_index as usize, v1_index as usize)) {
							Some(index) => *index as i32,
							None => -1,
						},
			};
		}
		-1
	}
	
	fn add_triangle(&mut self, mut v1: (i32, i32), mut v2: (i32, i32), mut v3: (i32, i32) ) {
		let ei = [self.has_edge(v1, v2), self.has_edge(v2, v3), self.has_edge(v3, v1)];
		
		if ei[0] >= 0 {
			let e_temp = self.get_edge(ei[0] as usize);
			v1 = self.get_vertex(e_temp.get_end()).pos();
			v2 = self.get_vertex(e_temp.get_start()).pos();
		} else if ei[1] >= 0 {
			let e_temp = self.get_edge(ei[1] as usize);
			v2 = self.get_vertex(e_temp.get_end()).pos();
			v3 = self.get_vertex(e_temp.get_start()).pos();
		} else if ei[2] >= 0 {
			let e_temp = self.get_edge(ei[2] as usize);
			v3 = self.get_vertex(e_temp.get_end()).pos();
			v1 = self.get_vertex(e_temp.get_start()).pos();
		}
		
		let vertecies = [self.get_or_add_vertex(v1.0, v1.1),
						 self.get_or_add_vertex(v2.0, v2.1),
						 self.get_or_add_vertex(v3.0, v3.1)];
						 
		let edges_len = self.edges.len();
		
		for index in 0..3 {
			let index_clamped = (index + 1) % 3;
			let next_edge_index = index_clamped + edges_len;
			let v_start = vertecies[index];
			let v_end = vertecies[index_clamped];
			let mut edge = Edge::new(v_start, v_end, next_edge_index);
			edge.set_opposite(ei[index]);
			if ei[index] >= 0 {
				self.get_edge_mut(ei[index] as usize).set_opposite((edges_len + index) as i32);
			}
			
			
			self.edges.push(edge);
			self.edge_map.insert((v_start, v_end), edges_len + index);
		}
		let tri = Triangle::new(edges_len);
		self.triangles.push(tri);
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
