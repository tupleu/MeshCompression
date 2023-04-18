use std::rc::{Rc, Weak};
use std::cell::RefCell;
use image::DynamicImage;
use std::collections::HashMap;
use rand::Rng;
use std::cmp;

const EPSILON: f32 = 0.0001_f32;

#[repr(C)]
#[derive(Copy, Clone, Debug, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub(crate) position: [f32; 3],
	pub(crate) color: [f32; 3],
	index: u32,
    anchor: u16, // 0: default, 1: x anchor, 2: y anchor, 3: xy anchor
	state: u16,
}

impl Vertex {
    fn new(position: [f32; 3], color: [f32; 3], index: u32, anchor: u16) -> Self {
        Vertex {
            position,
            color,
			index,
            anchor,
			state: 1,
        }
    }

    fn distance(a: Vertex, b: Vertex) -> f32 {
        let dx = f32::powf(a.position[0] + b.position[0], 2.0);
        let dy = f32::powf(a.position[1] + b.position[1], 2.0);
        let dz = f32::powf(a.position[2] + b.position[2], 2.0);
        f32::sqrt(dx + dy + dz)
    }
    fn midpoint(a: Vertex, b: Vertex) -> [f32; 3] {
        let mid_x = (a.position[0] + b.position[0]) / 2.0;
        let mid_y = (a.position[1] + b.position[1]) / 2.0;
        let mid_z = (a.position[2] + b.position[2]) / 2.0;
        [mid_x, mid_y, mid_z]
    }
    fn y_axis_aligned(a: Vertex, b: Vertex) -> bool {
        (a.position[0]-b.position[0]).abs() < EPSILON 
    }
    fn x_axis_aligned(a: Vertex, b: Vertex) -> bool {
        (a.position[1]-b.position[1]).abs() < EPSILON 
    }
    fn average_color(a: Vertex, b: Vertex) -> [f32; 3] {
        [(a.color[0] + b.color[0]) / 2.0, (a.color[1] + b.color[1]) / 2.0, (a.color[2] + b.color[2]) / 2.0]
    }
}

#[derive(Debug)]
pub struct Edge {
    vertex: Option<Rc<RefCell<Vertex>>>,
    opposite: Option<Rc<RefCell<Edge>>>,
    next: Option<Rc<RefCell<Edge>>>,
    triangle: Option<Rc<RefCell<Triangle>>>,
	history: Vec<Rc<RefCell<Vertex>>>,
	frozen: bool,
}

impl Edge {
    fn new() -> Rc<RefCell<Self>> {
        Rc::new(RefCell::new(Self {
            vertex: None,
            opposite: None,
            next: None,
            triangle: None,
			history: Vec::new(),
			frozen: false,
        }))
    }
    
    // fn length(self) -> Option<f32> {
    //     let a = self.vertex?.borrow().clone();
    //     let b = self.next?.clone().borrow().vertex?.clone().borrow();
    //     Some(Vertex::distance(a, *b))
    // }
}


fn distance(p1: [f32; 3], p2: [f32; 3]) -> f32 {
	let dx = f32::powf(p1[0] + p2[0], 2.0);
	let dy = f32::powf(p1[1] + p2[1], 2.0);
	let dz = f32::powf(p1[2] + p2[2], 2.0);
	f32::sqrt(dx + dy + dz)
}

#[derive(Debug)]
pub struct Triangle {
    edge: Option<Rc<RefCell<Edge>>>,
}

#[derive(Debug)]
pub struct Mesh {
    vertices: Vec<Rc<RefCell<Vertex>>>,
    edges: Vec<Rc<RefCell<Edge>>>,
    triangles: Vec<Rc<RefCell<Triangle>>>,
	vertex_edge_map: HashMap<(usize, usize), Rc<RefCell<Edge>>>,
    // edge_priority_queue,
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

	fn new_vertex(position: [f32; 3], color: [f32; 3], index: u32, anchor: u16) -> Rc<RefCell<Vertex>> {
		Rc::new(RefCell::new(Vertex::new(position, color, index, anchor)))
	}
	
	fn new_triangle(edge: Option<Rc<RefCell<Edge>>>) -> Rc<RefCell<Triangle>> {
		Rc::new(RefCell::new(Triangle { edge }))
	}
	
	//pub fn outline() -> Mesh {
		// same color, then make them closer together
		// specifics for edge collapse
		
		// float 16 to float 32
		
		
	//}
    pub fn length(edge: &Rc<RefCell<Edge>>) -> f32 {
      let v1 = *Mesh::vertex(edge).borrow();
      let v2 = *Mesh::vertex(&Mesh::next(&Mesh::next(edge))).borrow();
      Vertex::distance(v1, v2)
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
            let anchor = ((x == 0 || x == (width-1) as u32) as u16)*1 + ((y == 0 || y == (height-1) as u32) as u16)*2;
            println!("{},{},{}",x,y,anchor);
			mesh.vertices.push(Mesh::new_vertex([vx, vy*-1.0, 0.0], pixel.0, mesh.vertices.len() as u32, anchor));
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

    pub fn tri_count(&self) -> usize {
        self.triangles.len()
    }

	pub fn extract_vertices(&self) -> Vec<Vertex> {
		let mut vertices: Vec<Vertex> = Vec::new();
		for vertex in &self.vertices {
			let v: Vertex = *vertex.borrow();
			vertices.push(Vertex::new(v.position, v.color, 0, v.anchor));
		}
		// println!("Extracted verticies!");
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
		// println!("Extracted indicies!");
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
	
	pub fn triangle(edge: &Rc<RefCell<Edge>>) -> Rc<RefCell<Triangle>> {
		edge.borrow().triangle.clone().unwrap()
	}
	
	pub fn color_diff(edge: &Rc<RefCell<Edge>>) -> [f32; 3] {
		let c1 = Mesh::vertex(edge).borrow().color;
		let c2 = Mesh::vertex(&Mesh::next(edge)).borrow().color;
		[(c2[0] - c1[0]).abs(), (c2[1] - c1[1]).abs(), (c2[2] - c1[2].abs())]
	}
	
	pub fn get_neighborhood(&self, start_edge: &Rc<RefCell<Edge>>) -> Vec<Rc<RefCell<Edge>>> {
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
		//if edges.len() != 4 { println!("{}", edges.len()) }
		edges
	}// Loop <function> until <function>
	
	pub fn collapse_edge(&mut self) -> Result<(), String> {
        self.edges.sort_by(|a, b| Mesh::length(a).partial_cmp(&Mesh::length(b)).unwrap_or(std::cmp::Ordering::Equal));
        for edge in self.edges.iter() {
			// Check that the edge is not a boundary edge.
			// Todo: this is not always an invalid collapse (if both verticies are anchored in a shared direction), but will some modifications to deal with
			//if has_opposite {
			//    return Err("invalid collapse".to_string());
			//}
			
			if edge.borrow().frozen { return Err("invalid collapse".to_string()); } 
			
			
			let v1 = *Mesh::vertex(&edge).borrow();
			let v2 = *Mesh::vertex(&Mesh::next(&Mesh::next(&edge))).borrow();
			
			// last chance to invalidate
			// invalidate if there is not an opposite and the states are not equal or 0
			let has_opposite = Mesh::opposite(&edge).is_some();
			if !has_opposite && v1.state != 0 && v2.state != 0 && v1.state != v2.state {
				return Err("invalid collapse".to_string());
			}
			
			// Calculate the new position of the vertex.
            let a1 = v1.anchor;
            let a2 = v2.anchor;

            let p1 = v1.position;
            let p2 = v2.position;

            let position = match (a1, a2) {
                (0,0) => Vertex::midpoint(v1, v2),
                (3,0) => p1,
                (0,3) => p2,
                (1,0) | (0,1) | (1,1) if Vertex::y_axis_aligned(v1, v2) => Vertex::midpoint(v1, v2),
                (3,1) if Vertex::y_axis_aligned(v1, v2) => p1,
                (1,3) if Vertex::y_axis_aligned(v1, v2) => p2,
                (2,0) | (0,2) | (2,2) if Vertex::x_axis_aligned(v1, v2) => Vertex::midpoint(v1, v2),
                (3,2) if Vertex::x_axis_aligned(v1, v2) => p1,
                (2,3) if Vertex::x_axis_aligned(v1, v2) => p2,
                _ => continue,
            };
		
            let color = Vertex::average_color(v1, v2);
			
			let v_new = Mesh::new_vertex(position, color, self.vertices.len() as u32, a1 | a2);
			
			v_new.borrow_mut().state = if v1.state > v2.state { v1.state } else { v2.state } + 1;
			
			// Update all edges that use those verticies with v_new
			for e in self.get_neighborhood(&edge) {
				e.borrow_mut().vertex = Some(v_new.clone());
				e.borrow_mut().history.push(Mesh::vertex(&edge));
			}
			
			// do the same for the opposite, if it exists
			if has_opposite { 
				for e in self.get_neighborhood(&Mesh::opposite(&edge).unwrap()) {
					e.borrow_mut().vertex = Some(v_new.clone());
					e.borrow_mut().history.push(Mesh::vertex(&Mesh::next(&Mesh::next(&edge))));
				}
			}
			
			self.vertices.push(v_new);
			
			self.remove_triangle(Mesh::triangle(&edge)); 
			self.repair_opposites(&Mesh::next(&edge).borrow().opposite, &Mesh::next(&Mesh::next(&edge)).borrow().opposite);
			
			// do the same if we have an opposite
			// if has_opposite {
			// 	let opposite_edge = Mesh::opposite(&edge).unwrap();
			// 	self.remove_triangle(Mesh::triangle(&opposite_edge)); 
			// 	self.repair_opposites(&Mesh::next(&opposite_edge).borrow().opposite, &Mesh::next(&Mesh::next(&opposite_edge)).borrow().opposite);
			// }
		
			return Ok(());
		}
        Err("no valid collapses".to_string())
	}
	
	fn undo_edge_collapse(&mut self, ref_edge: &Rc<RefCell<Edge>>) {
		
		let edges = self.get_neighborhood(ref_edge);
		
		for edge in &edges {
			
			println!("{:?}", edge.borrow().history);
		}
		
	}
	
	pub fn undo_nearest_edge_collapse(&mut self, x: f32, y: f32) {
		let mut nearest_vertex = None;
		let point = [x, y, 0.0];
		let mut nearest = f32::INFINITY;
		for vertex in &self.vertices {
			if vertex.borrow().state < 2 { continue; }
			let d = distance(vertex.borrow().position, point);
			if d < nearest {
				nearest_vertex = Some(vertex);
				nearest = d;
			}
		}
		
		self.undo_edge_collapse_at(nearest_vertex.unwrap().clone());	
	}
	
	pub fn undo_edge_collapse_at(&mut self, at: Rc<RefCell<Vertex>>) {
		let mut found_edge = None;
		//let (v_remove, (v_revert1, v_revert2)) = self.history.remove(index as usize);
		for edge in &mut self.edges {
			if edge.borrow().frozen { continue; }
			if Rc::ptr_eq(&Mesh::vertex(edge), &at) {
				found_edge = Some(edge.clone());
				break;				
			}
		}
		self.undo_edge_collapse(&found_edge.unwrap());
	}
	
	fn repair_opposites(&mut self, e1: &Option<Rc<RefCell<Edge>>>, e2: &Option<Rc<RefCell<Edge>>>) {
		if e1.is_some() &&  e2.is_some() { 
			e1.clone().unwrap().borrow_mut().opposite = e2.clone(); 
			e2.clone().unwrap().borrow_mut().opposite = e1.clone();
			return;
		}
		if e1.is_none() &&  e2.is_some() { 
			e2.clone().unwrap().borrow_mut().opposite = None;
			return;
		}
		if e1.is_some() &&  e2.is_none() { 
			e1.clone().unwrap().borrow_mut().opposite = None; 
			return;
		}
		
	}
	
	pub fn remove_triangle(&mut self, triangle: Rc<RefCell<Triangle>>) -> bool {
		let edge_to_remove = Mesh::edge(&triangle);
		
		let mut index_to_remove: Option<usize> = None;
		for (i, triangle) in self.triangles.iter().enumerate() {
			if Rc::ptr_eq(&Mesh::edge(triangle), &edge_to_remove) {
				index_to_remove = Some(i);
				break;
			}
		}
		if index_to_remove.is_none() { return false; }
		
		// Freeze the edges
		edge_to_remove.borrow_mut().frozen = true;
		Mesh::next(&edge_to_remove).borrow_mut().frozen = true;
		Mesh::next(&Mesh::next(&edge_to_remove)).borrow_mut().frozen = true;
		
		self.triangles.remove(index_to_remove.unwrap());
		true
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
		
		e1.borrow_mut().triangle = self.triangles.last().cloned();
		e2.borrow_mut().triangle = self.triangles.last().cloned();
		e3.borrow_mut().triangle = self.triangles.last().cloned();
		
		self.vertex_edge_map.insert((indices[0], indices[1]), Rc::clone(&e1));
		self.vertex_edge_map.insert((indices[1], indices[2]), Rc::clone(&e2));
		self.vertex_edge_map.insert((indices[2], indices[0]), Rc::clone(&e3));
		
		self.edges.push(e1);
		self.edges.push(e2);
		self.edges.push(e3);
    }
}
