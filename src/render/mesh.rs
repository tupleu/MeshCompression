use std::rc::{Rc, Weak};
use std::cell::RefCell;
use image::DynamicImage;
use std::collections::HashMap;
use rand::Rng;

#[repr(C)]
#[derive(Copy, Clone, Debug, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub(crate) position: [f32; 3],
	pub(crate) color: [f32; 3],
	index: u32,
}

#[derive(Debug)]
pub struct Edge {
    vertex: Option<Rc<RefCell<Vertex>>>,
    opposite: Option<Rc<RefCell<Edge>>>,
    next: Option<Rc<RefCell<Edge>>>,
    triangle: Option<Rc<RefCell<Triangle>>>,
}

#[derive(Debug)]
pub struct Triangle {
    edge: Option<Rc<RefCell<Edge>>>,
}

impl Vertex {
    fn new(position: [f32; 3], color: [f32; 3], index: u32) -> Self {
        Vertex {
            position,
            color,
			index,
        }
    }
}

impl Edge {
    fn new() -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(Self {
            vertex: None,
            opposite: None,
            next: None,
            triangle: None,
        }))
    }
}

fn midpoint(p1: [f32; 3], p2: [f32; 3]) -> [f32; 3] {
	let mid_x = (p1[0] + p2[0]) / 2.0;
	let mid_y = (p1[1] + p2[1]) / 2.0;
	let mid_z = (p1[2] + p2[2]) / 2.0;
	[mid_x, mid_y, mid_z]
}

pub fn average_color(c1: &[f32; 3], c2: &[f32; 3]) -> [f32; 3] {
	[(c1[0] + c2[0]) / 2.0, (c1[1] + c2[1]) / 2.0, (c1[2] + c2[2]) / 2.0]
}

#[derive(Debug)]
pub struct Mesh {
    vertices: Vec<Rc<RefCell<Vertex>>>,
    edges: Vec<Rc<RefCell<Edge>>>,
    triangles: Vec<Rc<RefCell<Triangle>>>,
	vertex_edge_map: HashMap<(usize, usize), Rc<RefCell<Edge>>>,
}

impl Mesh {
    fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
			vertex_edge_map: HashMap::new(),
        }
    }

	fn new_vertex(position: [f32; 3], color: [f32; 3], index: u32) -> Rc<RefCell<Vertex>> {
		Rc::new(RefCell::new(Vertex::new(position, color, index)))
	}
	
	fn new_triangle(edge: Option<Rc<RefCell<Edge>>>) -> Rc<RefCell<Triangle>> {
		Rc::new(RefCell::new(Triangle { edge }))
	}

	pub fn from_image(dynamic_image: DynamicImage) -> Mesh {
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		
		// Create vertices
		for (x, y, pixel) in image.enumerate_pixels() {
			let vx = (x as f32 / max_dimension as f32) - 1.0;
			let vy = (y as f32 / max_dimension as f32) - 1.0;
			mesh.vertices.push(Mesh::new_vertex([vx, vy*-1.0, 0.0], pixel.0, mesh.vertices.len() as u32));
		}
		// Create triangles
		for y in 0..height-1 {
			for x in 0..width-1 {
				let i1 = x + (width * y);
				let i2 = x + (width * y) + 1;
				let i3 = x + (width * (y+1)) + 1;
				let i4 = x + (width * (y+1));
				mesh.add_triangle(&[i3, i2, i4]);
				mesh.add_triangle(&[i2, i1, i4]);
			}
		}
		mesh
	}

	pub fn extract_vertices(&self) -> Vec<Vertex> {
		let mut vertices: Vec<Vertex> = Vec::new();
		for vertex in &self.vertices {
			let v: Vertex = *vertex.borrow();
			vertices.push(Vertex::new(v.position, v.color, 0));
		}
		println!("Extracted verticies!");
		vertices
	}
	
	pub fn extract_indices(&self) -> Vec<u32> {
		let mut indices = Vec::new();
		for tri in &self.triangles {
			let mut current_edge = Rc::clone(&tri.borrow().edge.as_ref().unwrap());
			
			indices.push(Mesh::vertex(&current_edge).borrow().index);
			current_edge = Mesh::next(&current_edge);
			
			indices.push(Mesh::vertex(&current_edge).borrow().index);
			current_edge = Mesh::next(&current_edge);
			
			indices.push(Mesh::vertex(&current_edge).borrow().index);
		}
		println!("Extracted indicies!");
		indices
	}
	
	pub fn next(edge: &Rc<RefCell<Edge>>) -> Rc<RefCell<Edge>> {
		edge.borrow().next.clone().unwrap()
	}
	
	pub fn opposite(edge: &Rc<RefCell<Edge>>) -> Option<Rc<RefCell<Edge>>> {
		Some(edge.borrow().opposite.as_ref()?.clone())
	}
	
	pub fn vertex(edge: &Rc<RefCell<Edge>>) -> Rc<RefCell<Vertex>> {
		edge.borrow().vertex.clone().unwrap()
	}
	
	pub fn edge(triangle: &Rc<RefCell<Triangle>>) -> Rc<RefCell<Edge>> {
		triangle.borrow().edge.clone().unwrap()
	}
	
	pub fn get_neighboorhood(&self, start_edge: &Rc<RefCell<Edge>>) -> Vec<Rc<RefCell<Edge>>> {
		let mut edge = start_edge.clone();
		let mut prev_edge;
		let mut edges = Vec::new();
		loop {
			edges.push(edge.clone());
			prev_edge = edge.clone();
			let edge_opt = Mesh::opposite(&edge);
			if edge_opt.is_none() { break; }
			edge = Mesh::next(&edge_opt.unwrap());
			if Rc::ptr_eq(&edge, &start_edge) { break; }
		}
		edge = start_edge.clone();
		let mut edge_opt;
		loop {
			edge = Mesh::next(&Mesh::next(&edge));
			edge_opt = Mesh::opposite(&edge);
			if edge_opt.is_none() { break; }
			edge = edge_opt.unwrap();
			if Rc::ptr_eq(&edge, &prev_edge) { break; }
			if Rc::ptr_eq(&edge, &start_edge) { break; }
			edges.push(edge.clone());
			
		}
		
		edges
	}// Loop <function> until <function>
	
	pub fn collapse_edge(&mut self, edge: Rc<RefCell<Edge>>) {
		// Check that the edge is not a boundary edge.
		if Mesh::opposite(&edge).is_none() { return; }
		
		let opposite_edge = Mesh::opposite(&edge).unwrap()	;
		
		let v1 = Mesh::vertex(&edge);
		let v2 = Mesh::vertex(&opposite_edge);

		// Calculate the new position of the vertex.
		let position = midpoint(v1.borrow().position, v2.borrow().position);
		let color = average_color(&v1.borrow().color, &v2.borrow().color);
		let v_new = Mesh::new_vertex(position, color, self.vertices.len() as u32);
		
		// Update all edges that use those verticies with v_new
		for edge_ in self.get_neighboorhood(&edge) {
			edge_.borrow_mut().vertex = Some(v_new.clone());	
		}
		for edge_ in self.get_neighboorhood(&opposite_edge) {
			edge_.borrow_mut().vertex = Some(v_new.clone());
		}
		self.vertices.push(v_new);
	}
	
	pub fn get_random_edge(&self) -> Rc<RefCell<Edge>> {
		self.edges[rand::thread_rng().gen_range(0..self.edges.len() - 1)].clone()
	}

	fn find_vertex(&self, vertex: &Rc<RefCell<Vertex>>) -> usize {
		self.vertices.iter().position(|v| Rc::ptr_eq(v, vertex)).unwrap()
	}
	
	fn find_opposite_edge(&self, start: usize, end: usize) -> Option<Rc<RefCell<Edge>>> {
		self.vertex_edge_map.get(&(end, start)).cloned()
	}

    fn add_triangle(&mut self, indices: &[usize; 3]) {
		let v1 = self.vertices[indices[0]].clone();
		let v2 = self.vertices[indices[1]].clone();
		let v3 = self.vertices[indices[2]].clone();
		
		let e1 = Edge::new();
		let e2 = Edge::new();
		let e3 = Edge::new();
		
		e1.borrow_mut().vertex = Some(v1);
		e2.borrow_mut().vertex = Some(v2);
		e3.borrow_mut().vertex = Some(v3);
		
		e1.borrow_mut().next = Some(e2.clone());
		e2.borrow_mut().next = Some(e3.clone());
		e3.borrow_mut().next = Some(e1.clone());
		
		let e1_o = self.find_opposite_edge(indices[0], indices[1]);
		let e2_o = self.find_opposite_edge(indices[1], indices[2]);
		let e3_o = self.find_opposite_edge(indices[2], indices[0]);
		
		if e1_o.is_some() {
			e1.borrow_mut().opposite = Some(e1_o.clone().unwrap().clone());
			e1_o.unwrap().borrow_mut().opposite = Some(e1.clone());
		}
		if e2_o.is_some() {
			e2.borrow_mut().opposite = Some(e2_o.clone().unwrap().clone());
			e2_o.unwrap().borrow_mut().opposite = Some(e2.clone());
		}
		if e3_o.is_some() {
			e3.borrow_mut().opposite = Some(e3_o.clone().unwrap().clone());
			e3_o.unwrap().borrow_mut().opposite = Some(e3.clone());
		}
		
        self.triangles.push(Mesh::new_triangle(Some(e1.clone())));
		
		self.vertex_edge_map.insert((indices[0], indices[1]), Rc::clone(&e1));
		self.vertex_edge_map.insert((indices[1], indices[2]), Rc::clone(&e2));
		self.vertex_edge_map.insert((indices[2], indices[0]), Rc::clone(&e3));
		
		self.edges.push(e1);
		self.edges.push(e2);
		self.edges.push(e3);
    }
	
	
}