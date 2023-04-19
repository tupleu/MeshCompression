use std::rc::{Rc, Weak};
use std::cell::RefCell;
use image::DynamicImage;
use std::collections::HashMap;
use rand::Rng;
use std::cmp;

type Index = u32;
type Vec3f = [f32; 3];

const ANCHOR_NONE: u16 = 0;
const ANCHOR_X: u16 = 1;
const ANCHOR_Y: u16 = 2;
const ANCHOR_BOTH: u16 = 3;



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Helper Functions
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
fn average(p1: [f32; 3], p2: [f32; 3]) -> [f32; 3] {
	let mid_x = (p1[0] + p2[0]) / 2.0;
	let mid_y = (p1[1] + p2[1]) / 2.0;
	let mid_z = (p1[2] + p2[2]) / 2.0;
	[mid_x, mid_y, mid_z]
}



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Vertex
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub(crate) position: Vec3f,
	pub(crate) color: Vec3f,
	index: Index,
	anchor: u16,
	state: u16,
}

impl Vertex {
	pub fn new(position: Vec3f, color: Vec3f) -> Self {
		Self {
			position, 
			color, 
			index: 0,
			anchor: 0,
			state: 0,
		}
	}
}

#[derive(Debug, Clone)]
pub struct VertexPointer {
	vertex: Rc<RefCell<Vertex>>,
}

impl VertexPointer {
	// Constructors ///////////////////////////////////////////////////////////////////////////
    pub fn new(position: Vec3f, color: Vec3f, index: usize) -> Self {
        Self {  
			vertex: 
				Rc::new(RefCell::new(Vertex {
					position,
					color,
					index: index.try_into().unwrap(),
					anchor: ANCHOR_NONE,
					state: 1,
				}))
		}
    }
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn anchor(&self) -> u16 { self.vertex.borrow().anchor }
	pub fn state(&self) -> u16 { self.vertex.borrow().state }
	pub fn index(&self) -> Index { self.vertex.borrow().index }
	pub fn pos(&self) -> Vec3f { self.vertex.borrow().position }
	pub fn position(&self) -> Vec3f { self.pos() }
	pub fn color(&self) -> Vec3f { self.vertex.borrow().color }
	
	pub fn blend(&self, v: &VertexPointer) -> Vec3f {
		let c1 = self.color();
		let c2 = v.color();
		[(c1[0] + c2[0]) / 2.0, (c1[1] + c2[1]) / 2.0, (c1[2] + c2[2]) / 2.0]
	}
	
	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn set_color(&self, new_color: Vec3f) { self.vertex.borrow_mut().color = new_color; }
	pub fn set_state(&self, new_state: u16) { self.vertex.borrow_mut().state = new_state; }
	pub fn set_index(&self, new_index: Index) { self.vertex.borrow_mut().index = new_index; }
	pub fn set_anchor(&self, new_anchor: u16) { self.vertex.borrow_mut().anchor = new_anchor; }
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[derive(Debug)]
pub struct Edge {
    vertex: VertexPointer,
    opposite: Option<EdgePointer>,
    next: Option<EdgePointer>,
    triangle: Option<TrianglePointer>,
	history: Vec<VertexPointer>,
	index: Index, 
}

#[derive(Debug, Clone)]
pub struct EdgePointer {
	edge: Rc<RefCell<Edge>>,
}

impl EdgePointer {
	// Constructors ///////////////////////////////////////////////////////////////////////////
    fn new(vertex: &VertexPointer, index: usize) -> Self {
        Self {  
			edge: 
				Rc::new(RefCell::new(Edge {
					index: index.try_into().unwrap(),
					vertex: vertex.clone(),
					opposite: None,
					next: None,
					triangle: None,
					history: Vec::new(),
				}))
        }
    }
	// Getters ////////////////////////////////////////////////////////////////////////////////	
	pub fn next(&self) -> EdgePointer { self.edge.borrow().next.as_ref().unwrap().clone() }
	pub fn triangle(&self) -> TrianglePointer { self.edge.borrow().triangle.as_ref().unwrap().clone() }
	pub fn opposite(&self) -> Option<EdgePointer> { self.edge.borrow().opposite.clone() }
	pub fn has_opposite(&self) -> bool { self.edge.borrow().opposite.is_some() }
	pub fn vertex(&self) -> VertexPointer { self.edge.borrow().vertex.clone() }
	pub fn index(&self) -> Index { self.edge.borrow().index }
	pub fn neighboors(&self) -> Vec<EdgePointer> {
		let mut results = Vec::new();
		let mut current_edge = self.clone();
		loop {
			results.push(current_edge.clone());
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap().next();
			if current_edge.index() == self.index() { break; }

			println!("1? {:?}", current_edge.vertex().index());
		}
		if !current_edge.has_opposite() { return results; }
		current_edge = self.opposite().unwrap().next().next().clone();
		/*

		}

		loop {
			results.push(current_edge.clone());
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap().next();
			if current_edge.index() == self.index() { break; }
			println!("2? {:?}", current_edge.index());
		}
		*/
		results
	}
	
	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn set_index(&mut self, new_index: Index) { self.edge.borrow_mut().index = new_index; }
	pub fn set_vertex(&mut self, new_vertex: &VertexPointer) { self.edge.borrow_mut().vertex = new_vertex.clone(); }
	pub fn set_next(&mut self, new_next: &EdgePointer) { self.edge.borrow_mut().next = Some(new_next.clone()); }
	pub fn set_triangle(&mut self, new_triangle: &TrianglePointer) { self.edge.borrow_mut().triangle = Some(new_triangle.clone()); }
	pub fn set_opposite_some(&mut self, new_opposite: &EdgePointer) { self.set_opposite(&Some(new_opposite.clone())); }
	pub fn set_opposite(&mut self, new_opposite: &Option<EdgePointer>) { self.edge.borrow_mut().opposite = new_opposite.clone(); }
}



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[derive(Debug)]
pub struct Triangle {
    edge: EdgePointer,
	index: Index,
}

#[derive(Debug, Clone)]
pub struct TrianglePointer {
	triangle: Rc<RefCell<Triangle>>,
}

impl TrianglePointer {
	// Constructors ///////////////////////////////////////////////////////////////////////////
	fn new(edge: &EdgePointer, index: usize) -> Self {
        Self {  
			triangle: 
				Rc::new(RefCell::new(Triangle {
					index: index as Index,
					edge: edge.clone(),
				}))
        }
    }
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn edge(&self) -> EdgePointer { self.triangle.borrow().edge.clone() }
	pub fn index(&self) -> Index { self.triangle.borrow().index }
	
	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn set_index(&self, new_index: Index) { self.triangle.borrow_mut().index = new_index }
	

}

///////////////////////////////////////////////////////////////////////////////////////////////
// Mesh
///////////////////////////////////////////////////////////////////////////////////////////////

#[derive(Debug)]
pub struct Mesh {
    vertices: Vec<VertexPointer>,
    edges: Vec<EdgePointer>,
    triangles: Vec<TrianglePointer>,
	vertex_edge_map: HashMap<(Index, Index), EdgePointer>,
}

impl Mesh {
	// Constructors ///////////////////////////////////////////////////////////////////////////
	fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
			vertex_edge_map: HashMap::new(),
        }
    }
	
	pub fn from_image(dynamic_image: DynamicImage) -> Self { 
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		// Create vertices
		for (x, y, pixel) in image.enumerate_pixels() {
			
			let vx = (x as f32 / max_dimension as f32) - 1.0;
			let vy = (y as f32 / max_dimension as f32) - 1.0;
			let anchor = match (x as usize + 1, y as usize + 1) {
				(i, j) if (i == 1 || i == width) && (j == 1 || j == height) => ANCHOR_BOTH,
				(i, j) if (i == 1 || i == width) => ANCHOR_X,
				(i, j) if (j == 1 || j == height) => ANCHOR_Y,
				_ => ANCHOR_NONE,
			};
			let v = mesh.add_vertex([vx, vy*-1.0, 0.0], pixel.0);
			v.set_anchor(anchor);

		}
		println!("{}/{} vertices", mesh.vertices.len(), (width)*(height));
		// Create triangles
		for y in 0..height-1 {
			for x in 0..width-1 {
				let i1 = (x + (width * y)) as Index;
				let i2 = (x + (width * y) + 1) as Index;
				let i3 = (x + (width * (y+1)) + 1) as Index;
				let i4 = (x + (width * (y+1))) as Index;
				mesh.add_triangle(i1, i4, i2);
				mesh.add_triangle(i4, i3, i2);
			}
		}
		println!("{}/{} triangles", mesh.triangles.len(), (width-1)*(height-1)*2);
		println!("{}/{} edges", mesh.edges.len(), (width-1)*(height-1)*6);
		mesh
	}
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn triangle_count(&self) -> usize { self.triangles.len() }
	pub fn get_random_edge(&self) -> EdgePointer { self.edges[rand::thread_rng().gen_range(0..self.edges.len()-1)].clone() }
	pub fn extract_vertices(&self) -> Vec<Vertex> {
		let mut vertices = Vec::new();
		for vertex in &self.vertices {
			let v = vertex;
			vertices.push(Vertex::new(v.pos(), v.color()));
		}
		vertices
	}
	pub fn extract_indices(&self) -> Vec<Index> {
		let mut indices = Vec::new();
		for triangle in &self.triangles {
			indices.push(triangle.edge().vertex().index());			
			indices.push(triangle.edge().next().vertex().index());			
			indices.push(triangle.edge().next().next().vertex().index());			
		}
		indices
	}
	fn find_edge(&self, start: u32, end: u32) -> Option<EdgePointer> {
		self.vertex_edge_map.get(&(start, end)).cloned()
	}
	
	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn collapse_edge(&mut self, edge: EdgePointer) -> Result<(), String> {

		let (v1, v2) = (edge.vertex(), edge.next().next().vertex());
		let (a1, a2) = (v1.anchor(), v2.anchor());
		let (p1, p2) = (v1.position(), v2.position());
		
		let new_position = match (a1, a2) {
			(ANCHOR_BOTH, ANCHOR_BOTH) 	=>	return Err("Invalid collapse".to_string()),
			(ANCHOR_BOTH, ANCHOR_NONE) 	=>	p1,
			(ANCHOR_NONE, ANCHOR_BOTH) 	=>	p2,
			(_, ANCHOR_BOTH) |
				(ANCHOR_BOTH, _)		=>	return Err("Invalid collapse".to_string()),
			(ANCHOR_X, ANCHOR_Y) | 
				(ANCHOR_Y, ANCHOR_X) 	=>	return Err("Invalid collapse".to_string()),
			(ANCHOR_X, ANCHOR_NONE) 	=>	[p1[0], (p1[1] + p2[1]) / 2.0, 0.0],
			(ANCHOR_NONE, ANCHOR_X) 	=>	[p2[0], (p1[1] + p2[1]) / 2.0, 0.0],
			(ANCHOR_Y, ANCHOR_NONE) 	=>	[(p1[0] + p2[0]) / 2.0, p1[1], 0.0],
			(ANCHOR_NONE, ANCHOR_Y) 	=>	[(p1[0] + p2[0]) / 2.0, p2[1], 0.0],
			(ANCHOR_X, ANCHOR_X) 
				if p1[0] != p2[0] 		=>	return Err("Invalid collapse".to_string()),
			(ANCHOR_Y, ANCHOR_Y) 
				if p1[1] != p2[1]		=>	return Err("Invalid collapse".to_string()),
			_ 							=> 	average(p1, p2),
		};
		
		let new_color = average(v1.color(), v2.color());
		
		let v_new = self.add_vertex(new_position, new_color);

		v_new.set_state( if v1.state() > v2.state() { v1.state() + 1 } else { v2.state() + 1 } );
		
		// update all the edges with the new vertex
		for e in &mut edge.neighboors() {
			e.set_vertex(&v_new);
		}
		for e in &mut edge.next().neighboors() {
			e.set_vertex(&v_new);
		}
		
		// remove the triangles
		self.remove_triangle(&edge.triangle());
		if edge.has_opposite() {
			self.remove_triangle(&edge.opposite().unwrap().triangle());
		}
		
		// reconnect the opposites
		if edge.next().has_opposite() {
			edge.next().opposite().unwrap().set_opposite(&edge.next().next().opposite())
		}
		if edge.next().next().has_opposite() {
			edge.next().next().opposite().unwrap().set_opposite(&edge.next().opposite())
		}
		if edge.has_opposite() {
			let op_edge = edge.opposite().unwrap();
			if op_edge.next().has_opposite() {
				op_edge.next().opposite().unwrap().set_opposite(&op_edge.next().next().opposite())
			}
			if op_edge.next().next().has_opposite() {
				op_edge.next().next().opposite().unwrap().set_opposite(&op_edge.next().opposite())
			}
		}
		// Return that we finished properly
		Ok(())
	}
	
	fn add_vertex(&mut self, pos: Vec3f, color: Vec3f) -> VertexPointer { 
		self.vertices.push(VertexPointer::new(pos, color, self.vertices.len()));
		self.vertices.last().unwrap().clone()
	}
	fn add_edge(&mut self, start: &VertexPointer, end: &VertexPointer) -> EdgePointer {		self.edges.push(EdgePointer::new(start, self.edges.len()));
		self.vertex_edge_map.insert((start.index(), end.index()), self.edges.last().unwrap().clone());
		self.edges.last().unwrap().clone()
	}
	fn add_triangle(&mut self, i1: Index, i2: Index, i3: Index) {
		
		let v1 = self.vertices[i1 as usize].clone();
		let v2 = self.vertices[i2 as usize].clone();
		let v3 = self.vertices[i3 as usize].clone();
		
		let mut e1 = self.add_edge(&v1, &v2);
		let mut e2 = self.add_edge(&v2, &v3);
		let mut e3 = self.add_edge(&v3, &v1);
		
		e1.set_next(&e2);
		e2.set_next(&e3);
		e3.set_next(&e1);
		
		let e1_o = self.find_edge(i2, i1);
		let e2_o = self.find_edge(i3, i2);
		let e3_o = self.find_edge(i1, i3);
		
		if e1_o.is_some() {
			let mut opposite = e1_o.unwrap();
			e1.set_opposite_some(&opposite);
			opposite.set_opposite_some(&e1);
		}
		if e2_o.is_some() {
			let mut opposite = e2_o.unwrap();
			e2.set_opposite_some(&opposite);
			opposite.set_opposite_some(&e2);
		}
		if e3_o.is_some() {
			let mut opposite = e3_o.unwrap();
			e3.set_opposite_some(&opposite);
			opposite.set_opposite_some(&e3);
		} 
		
		self.triangles.push(TrianglePointer::new(&e1, self.triangles.len()));

		e1.set_triangle(self.triangles.last().unwrap());
		e2.set_triangle(self.triangles.last().unwrap());
		e3.set_triangle(self.triangles.last().unwrap());
	}
	
	fn remove_vertex(&mut self, vertex: &VertexPointer) -> bool {
		let index = vertex.index() as usize;
		println!("{} {} {}", vertex.index(), index, self.vertices.len());
		if index >= self.vertices.len().try_into().unwrap() { return false; }
		
		
		self.vertices.swap_remove(index);
		
		if index == self.vertices.len() { return false; }

		self.vertices[index].set_index(index as Index);
		
		true
	}
	fn remove_edge(&mut self, edge: &EdgePointer) -> bool {
		let index = edge.index() as usize;
		if index > self.edges.len().try_into().unwrap() { return false; }
		
		self.edges.swap_remove(index);
		self.edges[index].set_index(index as Index);
		
		true
	}
	fn remove_triangle(&mut self, triangle: &TrianglePointer) -> bool {
		let index = triangle.index() as usize;
		if index > self.triangles.len() { return false; }
		
		//self.remove_vertex(&triangle.edge().vertex());
		//self.remove_vertex(&triangle.edge().next().vertex());
		
		//self.remove_edge(&triangle.edge().next().next());
		//self.remove_edge(&triangle.edge().next());
		//self.remove_edge(&triangle.edge());
		
		//self.triangles.swap_remove(index);
		//self.triangles[index].set_index(index as Index);

		// Rc::strong_count()
		
		true
	}
}

/*
fn midpoint(p1: [f32; 3], p2: [f32; 3]) -> [f32; 3] {
	let mid_x = (p1[0] + p2[0]) / 2.0;
	let mid_y = (p1[1] + p2[1]) / 2.0;
	let mid_z = (p1[2] + p2[2]) / 2.0;
	[mid_x, mid_y, mid_z]
}

fn distance(p1: [f32; 3], p2: [f32; 3]) -> f32 {
	let dx = f32::powf(p1[0] + p2[0], 2.0);
	let dy = f32::powf(p1[1] + p2[1], 2.0);
	let dz = f32::powf(p1[2] + p2[2], 2.0);
	f32::sqrt(dx + dy + dz)
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

	pub fn from_image(dynamic_image: DynamicImage) -> Mesh {
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		//println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		
		// Create vertices
		let mut count = 0;
		for (x, y, pixel) in image.enumerate_pixels() {
			let vx = (x as f32 / max_dimension as f32) - 1.0;
			let vy = (y as f32 / max_dimension as f32) - 1.0;
            let corner = ((x == 0 || x == (width-1) as u32) && (y == 0 || y == (height-1) as u32)) as u16;
			mesh.vertices.push(Mesh::new_vertex([vx, vy*-1.0, 0.0], pixel.0, mesh.vertices.len() as u32, corner));
			//if corner != 0 { println!("{}", corner) }
			count += 1;
		}
		//println!("{} counted", count);
		//println!("{}/{} vertices", mesh.vertices.len(), (width)*(height));
		// Create triangles
		for y in 0..height-1 {
			for x in 0..width-1 {
				let i1 = x + (width * y);
				let i2 = x + (width * y) + 1;
				let i3 = x + (width * (y+1)) + 1;
				let i4 = x + (width * (y+1));
				mesh.add_triangle(&[i1, i4, i2]);
				mesh.add_triangle(&[i4, i3, i2]);
				
			}
		}
		//println!("{}/{} triangles", mesh.triangles.len(), (width-1)*(height-1)*2);
		//println!("{}/{} edges", mesh.edges.len(), (width-1)*(height-1)*6);
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
	
	fn is_clockwise(edge: &Rc<RefCell<Edge>>) -> bool {
		let v1 = Mesh::vertex(edge).borrow().position;
		let v2 = Mesh::vertex(&Mesh::next(edge)).borrow().position;
		let v3 = Mesh::vertex(&Mesh::next(&Mesh::next(edge))).borrow().position;
		
		let e1 = (v2[0] - v1[0])*(v2[0] + v1[0]);
		let e2 = (v3[0] - v2[0])*(v3[0] + v2[0]);
		let e3 = (v1[0] - v3[0])*(v1[0] + v3[0]);
		
		if e1+e2 < 0.0 || e2-e3 < 0.0{ return true; }
		
		false
	}
	
	pub fn extract_indices(&self) -> Vec<u32> {
		let mut count = 0;
		let mut indices = Vec::new();
		for tri in &self.triangles {
			let mut current_edge = Rc::clone(&tri.borrow().edge.as_ref().unwrap());
			count +=1;
			if current_edge.borrow().frozen {continue;}
			indices.push(Mesh::vertex(&current_edge).borrow().index);
			current_edge = Mesh::next(&current_edge);
			
			indices.push(Mesh::vertex(&current_edge).borrow().index);
			current_edge = Mesh::next(&current_edge);
			
			indices.push(Mesh::vertex(&current_edge).borrow().index);
			
		}
		//println!("Extracted {} triangles!", indices.len() / 3);
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
	
	pub fn get_neighboorhood(&self, start_edge: &Rc<RefCell<Edge>>) -> Vec<Rc<RefCell<Edge>>> {
		let mut edge = start_edge.clone();
		let mut prev_edge;
		let mut edges = Vec::new();
		loop {
			if !edge.borrow().frozen { 
				edges.push(edge.clone()); 
				//println!("	- {:?}", edge.borrow().vertex);
				////assert!!(Rc::ptr_eq(&Mesh::vertex(start_edge), &Mesh::vertex(&edge)));
			} else {
				//////assert!!(false);
				//println!("	~ {:?}", edge.borrow().vertex);
			}
			
			prev_edge = edge.clone();
			let edge_opt = Mesh::opposite(&edge);
			if edge_opt.is_none() { break; }
			edge = Mesh::next(&edge_opt.unwrap());
			if Rc::ptr_eq(&edge, &start_edge) { break; }
			continue;
			for e in &edges {
				////assert!!(!Rc::ptr_eq(&e, &edge));
			}
			
		}
		edge = start_edge.clone();
		/*
		let mut edge_opt;
		loop {
			
			edge = Mesh::next(&edge);
			edge_opt = Mesh::opposite(&edge);
			if edge_opt.is_none() { break; }
			edge = edge_opt.unwrap();
			if Rc::ptr_eq(&edge, &prev_edge) { break; }
			//if Rc::ptr_eq(&edge, &start_edge) { break; }
			if !edge.borrow().frozen { 
				edges.push(edge.clone());
				////assert!!(Rc::ptr_eq(&Mesh::vertex(start_edge), &Mesh::vertex(&edge)));
			}
			
		}*/
		//println!("length {}", edges.len());
		edges
	}// Loop <function> until <function>
	
	fn assert_frozen (&self, expected: usize) {
		return;
		let mut count = 0;
		for edge in &self.edges {
			if edge.borrow().frozen {continue;}
			let op = Mesh::opposite(&edge);
			if op.is_none() { continue; }
			if op.unwrap().borrow().frozen {count += 1;}
		}
		//println!("frozen count = {}/{}", count, expected);
		assert!(count <= expected);
	}
	
	pub fn collapse_edge(&mut self, edge: Rc<RefCell<Edge>>) -> Result<(), String> {

		//println!("collapse");
		self.assert_frozen(0);
		// Check that the edge is not a boundary edge.
		// Todo: this is not always an invalid collapse (if both verticies are anchored in a shared direction), but will some modifications to deal with
		//if has_opposite {
        //    return Err("invalid collapse".to_string());
        //}
		if Mesh::opposite(&edge).is_none() {
            return Err("invalid collapse".to_string());
        }
	
		
		if edge.borrow().frozen { return Err("invalid collapse".to_string()); } 
		
		let next_edge = &Mesh::next(&edge);
		
		let v1 = Mesh::vertex(&edge);
		let v2 = Mesh::vertex(&next_edge);
		//let v3 = Mesh::vertex(&Mesh::opposite(&edge).unwrap());
		//println!("{:?}", v2);
		//println!("{:?}", v3);
		
		if v1.borrow().anchor == 1 && v2.borrow().anchor == 1 {
            return Err("invalid collapse".to_string());
        }
		
		// last chance to invalidate
		// invalidate if there is not an opposite and the states are not equal or 0
		let has_opposite = Mesh::opposite(&edge).is_some();
		if !has_opposite && v1.borrow().state != 0 && v2.borrow().state != 0 && v1.borrow().state != v2.borrow().state {
			return Err("invalid collapse".to_string());
		}
		
		// Calculate the new position of the vertex.
        // println!("{:?}, {:?}", v1.borrow().anchor, v2.borrow().anchor);
		let position = match (v1.borrow().anchor, v2.borrow().anchor) {
		    (0,0) => midpoint(v1.borrow().position, v2.borrow().position),
		    (1,0) => v1.borrow().position,
		    (0,1) => v2.borrow().position,
            _ => return Err("invalid collapse".to_string()),
		};
		
		let color = average_color(&v1.borrow().color, &v2.borrow().color);
        
		let v_new = Mesh::new_vertex(position, color, self.vertices.len() as u32, v1.borrow().anchor | v2.borrow().anchor);
		
		v_new.borrow_mut().state = if v1.borrow().state > v2.borrow().state { v1.borrow().state } else { v2.borrow().state } + 1;
		
		
		self.remove_triangle(&edge); 
		if has_opposite {
			////assert!!(!Mesh::opposite(&edge).unwrap().borrow().frozen);
			////assert!!(!Rc::ptr_eq(&edge, &Mesh::opposite(&edge).unwrap()));
			////assert!!(!Rc::ptr_eq(&Mesh::vertex(&edge), &Mesh::vertex(&Mesh::opposite(&edge).unwrap())));
			self.assert_frozen(3);
			self.remove_triangle(&Mesh::opposite(&edge).unwrap()); 
			self.assert_frozen(4);
		}
		
		// Update all edges that use those verticies with v_new
		for e in self.get_neighboorhood(&edge) {
			
			//println!("{:?}", e.borrow().vertex);
			e.borrow_mut().vertex = Some(v_new.clone());
			e.borrow_mut().history.push(v1.clone());
			//println!("hi");
		}
		
		// do the same for the opposite, if it exists
		for e in self.get_neighboorhood(&next_edge) {
			//println!("{:?}", e.borrow().vertex);
			e.borrow_mut().vertex = Some(v_new.clone());
			e.borrow_mut().history.push(v2.clone());
			
			//println!("hi2");
		}
		
		self.vertices.push(v_new);
		
		
		
		
		
		
		self.repair_opposites(&Mesh::opposite(&Mesh::next(&edge)), &Mesh::opposite(&Mesh::next(&Mesh::next(&edge))));
		
		if has_opposite {
			self.assert_frozen(2);
			self.repair_opposites(&Mesh::opposite(&Mesh::next(&Mesh::opposite(&edge).unwrap())), &Mesh::opposite(&Mesh::next(&Mesh::next(&Mesh::opposite(&edge).unwrap()))));
		}
		
		self.assert_frozen(0);
		
		
        Ok(())
	}
	
	fn undo_edge_collapse(&mut self, ref_edge: &Rc<RefCell<Edge>>) {
		
		let edges = self.get_neighboorhood(ref_edge);
		//println!("-----");
		for edge in &edges {
			let v = edge.borrow_mut().history.pop();
			edge.borrow_mut().vertex = v;
			if edge.borrow().frozen {
				//println!("frozen");
			}
			//println!("{:?}", edge.borrow().history);
		}
		
	}
	
	pub fn undo_nearest_edge_collapse(&mut self, x: f32, y: f32) -> bool {
		let mut nearest_vertex = Vec::new();
		let point = [x, y, 0.0];
		for vertex in &mut self.vertices {
			if vertex.borrow().state < 2 { continue; }
			let d = distance(vertex.borrow().position, point);
			nearest_vertex.push((vertex.clone(), d));
		}
		nearest_vertex.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
		
		loop {
			if nearest_vertex.is_empty() { break; }
			if self.undo_edge_collapse_at(&nearest_vertex.pop().unwrap().0) {
				return true;
			}	
		}
		
		false
	}
	
	pub fn undo_edge_collapse_at(&mut self, at: &Rc<RefCell<Vertex>>) -> bool {
		let mut found_edge = None;
		//let (v_remove, (v_revert1, v_revert2)) = self.history.remove(index as usize);
		for edge in &mut self.edges {
			if edge.borrow().frozen { continue; }
			if Rc::ptr_eq(&Mesh::vertex(edge), at) {
				found_edge = Some(edge.clone());
				break;				
			}
		}
		if found_edge.is_none() { return false; }
		self.undo_edge_collapse(&found_edge.unwrap());
		true
	}
	
	fn repair_opposites(&mut self, e1: &Option<Rc<RefCell<Edge>>>, e2: &Option<Rc<RefCell<Edge>>>) {
		if e1.is_some() &&  e2.is_some() { 
			e1.clone().unwrap().borrow_mut().opposite = e2.clone(); 
			e2.clone().unwrap().borrow_mut().opposite = e1.clone();
			////assert!!(!Rc::ptr_eq(&e1.clone().unwrap(), &e2.clone().unwrap()));
			////assert!!(Rc::ptr_eq(&Mesh::opposite(&e1.clone().unwrap()).unwrap(), &e2.clone().unwrap()));
			////assert!!(Rc::ptr_eq(&Mesh::opposite(&e2.clone().unwrap()).unwrap(), &e1.clone().unwrap()));
			////assert!!(!e1.clone().unwrap().borrow().frozen);
			////assert!!(!e2.clone().unwrap().borrow().frozen);
			//println!("??");
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
	
	pub fn remove_triangle(&mut self, edge_to_remove: &Rc<RefCell<Edge>>) -> bool {
		//////assert!!(Rc::ptr_eq(&triangle, &Mesh::triangle(&edge_to_remove)));
		/*let mut index_to_remove: Option<usize> = None;
		for (i, t) in self.triangles.iter().enumerate() {
			if Rc::ptr_eq(&t, &Mesh::triangle(&edge_to_remove)) {
				index_to_remove = Some(i);
				break;
			}
		}
		if index_to_remove.is_none() { return false; }*/
		
		// Freeze the edges
		let edge_next = Mesh::next(&edge_to_remove);
		let edge_next_next = Mesh::next(&Mesh::next(&edge_to_remove));
		edge_to_remove.borrow_mut().frozen = true;
		edge_next.borrow_mut().frozen = true;
		edge_next_next.borrow_mut().frozen = true;
		////assert!!(Rc::ptr_eq(&Mesh::triangle(edge_to_remove), &Mesh::triangle(&edge_next)));
		////assert!!(Rc::ptr_eq(&Mesh::triangle(edge_to_remove), &Mesh::triangle(&edge_next_next)));
		//println!("~~~~~~~~~~~~~~~~~~~~~~~~~{:?}", Mesh::vertex(edge_to_remove));
		//self.triangles.remove(index_to_remove.unwrap());
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
		
        self.triangles.push(Mesh::new_triangle(Some(e3.clone())));
		
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
*/