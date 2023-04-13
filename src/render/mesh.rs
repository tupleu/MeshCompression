use std::collections::HashMap;
use image::DynamicImage;

pub mod vertex;
mod edge;
//mod triangle;
use self::vertex::{Point, Vertex, VertexIndex};
use self::edge::{Edge, EdgeIndex};
//use self::triangle::Triangle;



pub struct Triangle {
	v_start: usize,
	edge: EdgeIndex,
}
impl Triangle { 
	pub fn new(v_start: usize, edge: EdgeIndex) -> Triangle { Triangle { v_start: v_start, edge: edge } }
	pub fn edge(&self) -> EdgeIndex { self.edge }
}


pub struct Mesh {
	vertices: Vec<Vertex>,
	edges: Vec<Edge>,
	triangles: Vec<Triangle>,
	triangle_verticies: Vec<VertexIndex>,
	vertex_map: HashMap<Point, VertexIndex>,
	edge_map: HashMap<(VertexIndex, VertexIndex), EdgeIndex>,
}

impl Mesh {
	pub fn new() -> Mesh { Mesh { vertices: vec![], edges: vec![], triangles: vec![], triangle_verticies: vec![], vertex_map: HashMap::new(), edge_map: HashMap::new() } }
	
	pub fn from_image(dimg: DynamicImage) -> Mesh {
		let width = dimg.width() as usize;
		let height = dimg.height() as usize;
		let img = dimg.to_rgb32f();
		
		let mut mesh = Mesh::new();
		
		let mut vs = vec![vec![0; height+1]; width + 1];
		// make the verticies
		let max_d = if width > height { width } else { height } as f32;
		for (i, j, pixel) in img.enumerate_pixels() {
			vs[i as usize][j as usize] = mesh.vertices.len();
			mesh.m_get_or_add_vertex(2.0*(i as f32)/max_d-1.0, 1.0-2.0*(j as f32)/max_d, Some(pixel.0));
		}
		
		// create the edges and faces
		for (i, j, _) in img.enumerate_pixels() {
			let x = i as f32;
			let y = j as f32;
			mesh.m_add_triangle((2.0*(x as f32)/max_d-1.0, 1.0-2.0*(y as f32)/max_d),
								(2.0*(x as f32)/max_d-1.0, 1.0-2.0*(y+1 as f32)/max_d),
								(2.0*(x+1 as f32)/max_d-1.0, 1.0-2.0*(y as f32)/max_d));
			mesh.m_add_triangle((2.0*(x as f32)/max_d-1.0, 1.0-2.0*(y+1 as f32)/max_d),
								(2.0*(x+1 as f32)/max_d-1.0, 1.0-2.0*(y+1 as f32)/max_d),
								(2.0*(x+1 as f32)/max_d-1.0, 1.0-2.0*(y as f32)/max_d));
		}
		// Return the mesh		
		mesh
	}
	
	fn m_get_or_add_vertex(&mut self, x: f32, y: f32, color: Option<[f32; 3]>) -> VertexIndex {
		self.get_or_add_vertex(Point::new(x, y), color)
	}
	
	fn m_add_triangle(&mut self, v1: (f32, f32), v2: (f32, f32), v3: (f32, f32)) { 
		self.add_triangle(Point::new(v1.0, v1.1), Point::new(v2.0, v2.1), Point::new(v3.0, v3.1))
	}
	
	pub fn vertices(&self) -> &Vec<Vertex> {
		return &self.vertices;
	}
	
	pub fn indices(&self) -> Vec<u32> {
		self.triangle_verticies.iter().map(|t| t.index).collect()
	}
	
	fn get_or_add_vertex(&mut self, p: Point, color: Option<[f32; 3]>) -> VertexIndex {
		if let Some(index) = self.vertex_map.get(&p) {
			return *index;
		}
		let index = VertexIndex::new(self.vertices.len() as u32);
		self.vertices.push( Vertex::new(p.x, p.y, color.unwrap_or([0.0, 0.0, 0.0])) );
		self.vertex_map.insert(p, index);
		index
	}
	
	pub fn get_vertex(&self, index: VertexIndex) -> &Vertex {
		&self.vertices[index.index as usize]
	}
	
	pub fn get_edge(&self, index: EdgeIndex) -> &Edge {
		&self.edges[index.index as usize]
	}
	
	pub fn get_vertex_mut(&mut self, index: VertexIndex) -> &mut Vertex {
		&mut self.vertices[index.index as usize]
	}
	
	fn get_edge_mut(&mut self, index: EdgeIndex) -> &mut Edge {
		&mut self.edges[index.index as usize]
	}
	
	pub fn find_vertex(&self, v: Point) -> Option<VertexIndex> {
		self.vertex_map.get(&v).copied()
	}
	
	pub fn find_edge(&self, v_start: Point, v_end: Point) -> Option<EdgeIndex> {
		let v1_index = self.find_vertex(v_start);
		let v2_index = self.find_vertex(v_end);
		if v1_index.is_none() || v2_index.is_none() { return None; }
		
		let result = self.edge_map.get(&(v1_index.unwrap(), v2_index.unwrap()));
		if result.is_some() { return result.copied(); }
		
		self.edge_map.get(&(v2_index.unwrap(), v1_index.unwrap())).copied()
	}
	
	pub fn add_triangle(&mut self, mut v1: Point, mut v2: Point, mut v3: Point ) {
		let ei = [self.find_edge(v1, v2), self.find_edge(v2, v3), self.find_edge(v3, v1)];
		
		if ei[0].is_some() {
			let e_temp = self.get_edge(ei[0].unwrap());
			v1 = self.get_vertex(e_temp.end()).pos();
			v2 = self.get_vertex(e_temp.start()).pos();
		} else if ei[1].is_some() {
			let e_temp = self.get_edge(ei[1].unwrap());
			v2 = self.get_vertex(e_temp.end()).pos();
			v3 = self.get_vertex(e_temp.start()).pos();
		} else if ei[2].is_some() {
			let e_temp = self.get_edge(ei[2].unwrap());
			v3 = self.get_vertex(e_temp.end()).pos();
			v1 = self.get_vertex(e_temp.start()).pos();
		}
		
		let vertices = [self.get_or_add_vertex(v1, None),
						 self.get_or_add_vertex(v2, None),
						 self.get_or_add_vertex(v3, None)];
						 
		let edges_len = self.edges.len();
		
		for index in 0..3 {
			let index_clamped = (index + 1) % 3;
			let next_edge_index = index_clamped + edges_len;
			let v_start = vertices[index];
			let v_end = vertices[index_clamped];
			let mut edge = Edge::new(v_start, v_end, EdgeIndex::new(next_edge_index as u32));
			edge.set_opposite(ei[index]);
			if ei[index].is_some() {
				self.get_edge_mut(ei[index].unwrap()).set_opposite(Some(EdgeIndex::new((edges_len + index) as u32)));
			}
			
			
			self.edges.push(edge);
			self.edge_map.insert((v_start, v_end), EdgeIndex::new((edges_len + index) as u32));
		}
		self.triangles.push(Triangle::new(self.triangle_verticies.len(), EdgeIndex::new(self.edges.len() as u32 - 1)));
		self.triangle_verticies.push(vertices[0]);
		self.triangle_verticies.push(vertices[1]);
		self.triangle_verticies.push(vertices[2]);
		
	}

	pub fn remove_triangle(&mut self, triangle: Triangle) {}
	
	pub fn set_vertex_color(&mut self, vertex: VertexIndex, color: [f32; 3]) {}
	
	pub fn next(&self, edge: EdgeIndex) -> EdgeIndex { EdgeIndex::new(0) }
	
	pub fn opposite(&self, edge: EdgeIndex) -> Option<EdgeIndex> { None }
	
	pub fn start(&self, edge: EdgeIndex) -> Option<&Vertex> { None }
	
	pub fn length(&self, edge: EdgeIndex) -> f32 { 0.0 }

	/*
	fn add_triangle_to_edge(&mut self, edge_index: usize, vertex: (i32, i32)) {
		let vi1 = self.get_edge(edge_index).end();
		let vi2 = self.get_edge(edge_index).start();
		
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


