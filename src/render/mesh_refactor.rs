use std::rc::{Rc, Weak};
use std::cell::RefCell;
use image::DynamicImage;
use std::collections::{HashMap, HashSet, BTreeSet, BinaryHeap};
use rand::Rng;
use rand::seq::SliceRandom;
use std::cmp::{self, Reverse};

type Index = u32;
type Vec3f = [f32; 3];

const ANCHOR_NONE: u16 = 0;
const ANCHOR_X: u16 = 1;
const ANCHOR_Y: u16 = 2;
const ANCHOR_BOTH: u16 = 3;

const EPSILON: f32 = 0.001_f32;



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

fn distance(p1: [f32; 3], p2: [f32; 3]) -> f32 {
	let dx = f32::powf(p1[0] - p2[0], 2.0);
	let dy = f32::powf(p1[1] - p2[1], 2.0);
	let dz = f32::powf(p1[2] - p2[2], 2.0);
	//println!("{:?} {:?}, {} {}", p1, p2, dx, dy);
	f32::sqrt(dx + dy + dz)
}

// returns (m, b)
fn line_equation(edge: &EdgePointer) -> (f32, f32) {
	let p1 = edge.vertex().pos();
	let p2 = edge.next().vertex().pos();
	let m = (p2[1] - p1[1]) / (p2[0] - p1[0]);
	(m, p1[1] - (m * p1[0]))
}

fn in_range(i: f32, r: (f32, f32)) -> bool {
	let r1 = if r.0 > r.1 { r.1 } else { r.0 };
	let r2 = if r.0 > r.1 { r.0 } else { r.1 };
	if i > r1 + f32::EPSILON && i < r2 - f32::EPSILON { return true; }
	false
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Vertex
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[repr(C)]
#[derive(Copy, Clone, Default, Debug, PartialEq, bytemuck::Pod, bytemuck::Zeroable)]
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
			state: 1,
		}
	}
}

#[derive(Debug, Clone, PartialEq)]
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
	fn distance_to(&self, v: &VertexPointer) -> f32 { distance(self.pos(), v.pos()) }
	
	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn set_color(&self, new_color: Vec3f) { self.vertex.borrow_mut().color = new_color; }
	pub fn set_position(&self, new_position: Vec3f) { self.vertex.borrow_mut().position = new_position; }
	pub fn set_state(&self, new_state: u16) { self.vertex.borrow_mut().state = new_state; }
	pub fn set_index(&self, new_index: Index) { self.vertex.borrow_mut().index = new_index; }
	pub fn set_anchor(&self, new_anchor: u16) { self.vertex.borrow_mut().anchor = new_anchor; }
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[derive(Debug, PartialEq)]
pub struct Edge {
    vertex: VertexPointer,
    opposite: Option<EdgePointer>,
    next: Option<EdgePointer>,
    triangle: Option<TrianglePointer>,
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
					vertex: vertex.clone(),
					opposite: None,
					next: None,
					triangle: None,
					index: index.try_into().unwrap(),
				}))
        }
    }
	// Getters ////////////////////////////////////////////////////////////////////////////////	
	pub fn next(&self) -> EdgePointer { self.edge.borrow().next.as_ref().unwrap().clone() }
	pub fn triangle(&self) -> TrianglePointer { self.edge.borrow().triangle.as_ref().unwrap().clone() }
	pub fn opposite(&self) -> Option<EdgePointer> { self.edge.borrow().opposite.clone() }
	pub fn length(&self) -> f32 { self.vertex().distance_to(&self.next().vertex()) }
	pub fn has_opposite(&self) -> bool { self.edge.borrow().opposite.is_some() }
	pub fn vertex(&self) -> VertexPointer { self.edge.borrow().vertex.clone() }
	pub fn index(&self) -> Index { self.edge.borrow().index }
	pub fn intersects(&self, edge: &EdgePointer) -> bool {
		
		let r1 = (self.vertex().pos()[0], self.next().vertex().pos()[0]);
		let r2 = (edge.vertex().pos()[0], edge.next().vertex().pos()[0]);
		let (m1, b1) = line_equation(self);
		let (m2, b2) = line_equation(edge);
		if f32::abs(m2 - m1) < f32::EPSILON { return false; }
		
		let x = (b2 - b1) / (m1 - m2);
		let y = (x * m1) + b1;
		if in_range(x, r1) && in_range(x, r2) {
			return true;			
		}
		false 
	}
	pub fn neighbors(&self) -> Vec<EdgePointer> {
		let self_vi = self.vertex().index();
		let mut results = Vec::new();
		let mut current_edge = self.clone();
		loop {
			results.push(current_edge.clone());
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap().next();
			if current_edge.index() == self.index() { break; }
			assert_eq!(current_edge.vertex().index(), self_vi);
		}
		current_edge = self.next().next();
		if !current_edge.has_opposite() { return results; }
		current_edge = current_edge.opposite().unwrap();
		loop {
			results.push(current_edge.clone());
			current_edge = current_edge.next().next();
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap();
			if current_edge.index() == self.index() { break; }
			assert_eq!(current_edge.vertex().index(), self_vi);
		}
		results
	}
    pub fn energy(&self) -> f32 {
        let length: f32 = self.length();
        let color_diff: f32 = distance(self.vertex().color(), self.next().vertex().color());
        length + color_diff/(length+1.0)
    }
	

	// Setters ////////////////////////////////////////////////////////////////////////////////
	pub fn set_index(&mut self, new_index: Index) { self.edge.borrow_mut().index = new_index; }
	pub fn set_vertex(&mut self, new_vertex: &VertexPointer) { self.edge.borrow_mut().vertex = new_vertex.clone(); }
	pub fn set_next(&mut self, new_next: &EdgePointer) { self.edge.borrow_mut().next = Some(new_next.clone()); }
	pub fn set_triangle(&mut self, new_triangle: &TrianglePointer) { self.edge.borrow_mut().triangle = Some(new_triangle.clone()); }
	pub fn set_opposite_some(&mut self, new_opposite: &EdgePointer) { self.set_opposite(&Some(new_opposite.clone())); }
	pub fn set_opposite(&mut self, new_opposite: &Option<EdgePointer>) { self.edge.borrow_mut().opposite = new_opposite.clone(); }
}

impl Ord for EdgePointer {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
       self.energy().total_cmp(&other.energy()) 
    }
}
impl PartialOrd for EdgePointer {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for EdgePointer {
    fn eq(&self, other: &Self) -> bool {
        self.energy() == other.energy()
    }
}
impl Eq for EdgePointer {}




///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
#[derive(Debug, PartialEq)]
pub struct Triangle {
    edge: EdgePointer,
	index: Index,
	centroid: (f32, f32),
	radius: f32,
}

#[derive(Debug, Clone, PartialEq)]
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
					centroid: (0.0, 0.0),
					radius: 0.0,
				}))
        }
    }
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn edge(&self) -> EdgePointer { self.triangle.borrow().edge.clone() }
	pub fn index(&self) -> Index { self.triangle.borrow().index }
	pub fn edges(&self) -> [EdgePointer; 3] {
		[self.edge(), self.edge().next(), self.edge().next().next()]
	}
	pub fn area(&self) -> f32 {
		let p1 = self.edge().vertex().pos();
		let p2 = self.edge().next().vertex().pos();
		let p3 = self.edge().next().next().vertex().pos();
		0.5 * (
			(p1[0] * (p2[1] - p3[1])) +
			(p2[0] * (p3[1] - p1[1])) +
			(p3[0] * (p1[1] - p2[1]))
		)
	}
	pub fn circumscribed_radius(&self) -> f32 {
		self.triangle.borrow().radius
	}
	pub fn centroid(&self) -> Vec3f {
		let c = self.triangle.borrow().centroid;
		[c.0, c.1, 0.0]
	}
	pub fn update(&mut self) {
		self.update_centroid();
		self.update_circumscribed_radius();
	}
	pub fn update_centroid(&mut self) {
		let v1 = self.triangle.borrow().edge.vertex().pos();
		let v2 = self.triangle.borrow().edge.next().vertex().pos();
		let v3 = self.triangle.borrow().edge.next().next().vertex().pos();
		let sum_x = v1[0] + v2[0] + v3[0];
		let sum_y = v1[1] + v2[1] + v3[1];
		self.triangle.borrow_mut().centroid = (sum_x / 3.0, sum_y / 3.0);
	}
	pub fn update_circumscribed_radius(&mut self) {
		
		self.triangle.borrow_mut().radius = (self.edge().length() * self.edge().next().length() * self.edge().next().next().length()) / (4.0 * self.area());
		//println!("area = {}", self.edge().length());
	}
	pub fn overlaps(&self, triangle: &TrianglePointer) -> bool {
		// do a quick check for circular overlap
		
		let d = distance(self.centroid(), triangle.centroid());
		if d > self.circumscribed_radius() + triangle.circumscribed_radius() {
			//println!("d = {}", d);
			//println!("r = {}", self.circumscribed_radius());
			return false;
		}
		
		// Check each of the edges
		for edge1 in self.edges() {
            for edge2 in triangle.edges() {
                if edge1.intersects(&edge2) {
					//println!("uhh");
                    return true;
                }
            }
        }
		false
	}
	
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
	energy_set: BTreeSet<EdgePointer>,
	vertex_edge_map: HashMap<(Index, Index), EdgePointer>,
	history: Vec<(Index, Index, Index, (f32, f32), (f32, f32))>, // vl, vr, vs, (dx dy), (dx, dy)
}

impl Mesh {
	// Constructors ///////////////////////////////////////////////////////////////////////////
	fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
			energy_set: BTreeSet::new(),
			vertex_edge_map: HashMap::new(),
			history: Vec::new(), 
        }
    }
	
	pub fn from_image1(dynamic_image: DynamicImage) -> Self { 
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		// Create vertices
		for (x, y, pixel) in image.enumerate_pixels() {
			let vx = (x as f32) / max_dimension - 1.0;
			let vy = -1_f32*((y as f32) / max_dimension - 1.0);
			let anchor = match (x as usize + 1, y as usize + 1) {
				(i, j) if (i == 1 || i == width) && (j == 1 || j == height) => ANCHOR_BOTH,
				(i, j) if (i == 1 || i == width) => ANCHOR_X,
				(i, j) if (j == 1 || j == height) => ANCHOR_Y,
				_ => ANCHOR_NONE,
			};
			let v = mesh.add_vertex([vx, vy, 0.0], pixel.0);
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

	pub fn from_image2(dynamic_image: DynamicImage) -> Self { 
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		// Create vertices
		for (x, y, pixel) in image.enumerate_pixels() {
			let vx = (x as f32) / max_dimension - 1.0;
			let vy = -1_f32*((y as f32) / max_dimension - 1.0);
			let anchor = match (x as usize + 1, y as usize + 1) {
				(i, j) if (i == 1 || i == width) && (j == 1 || j == height) => ANCHOR_BOTH,
				(i, j) if (i == 1 || i == width) => ANCHOR_X,
				(i, j) if (j == 1 || j == height) => ANCHOR_Y,
				_ => ANCHOR_NONE,
			};
			let v = mesh.add_vertex([vx, vy, 0.0], pixel.0);
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
                if (x+y) & 1 == 1 {
                    mesh.add_triangle(i1, i4, i2);
                    mesh.add_triangle(i4, i3, i2);
                }
                else {
                    mesh.add_triangle(i1, i3, i2);
                    mesh.add_triangle(i4, i3, i1);
                }
			}
		}
		println!("{}/{} triangles", mesh.triangles.len(), (width-1)*(height-1)*2);
		println!("{}/{} edges", mesh.edges.len(), (width-1)*(height-1)*6);
		mesh
	}
	pub fn from_image3(dynamic_image: DynamicImage) -> Self { 
        todo!("not finished");
		// let width = dynamic_image.width() as usize;
		// let height = dynamic_image.height() as usize;
		// println!("{}x{}", width, height);
		// let image = dynamic_image.to_rgb32f();
		//
		// let mut mesh = Mesh::new();
		//
		// let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;
		// // Create vertices
		// for (x, y, pixel) in image.enumerate_pixels() {
		// 	let vx = (x as f32) / max_dimension - 1.0;
		// 	let vy = -1_f32*((y as f32) / max_dimension - 1.0);
		// 	let anchor = match (x as usize + 1, y as usize + 1) {
		// 		(i, j) if (i == 1 || i == width) && (j == 1 || j == height) => ANCHOR_BOTH,
		// 		(i, j) if (i == 1 || i == width) => ANCHOR_X,
		// 		(i, j) if (j == 1 || j == height) => ANCHOR_Y,
		// 		_ => ANCHOR_NONE,
		// 	};
		// 	let v = mesh.add_vertex([vx, vy, 0.0], pixel.0);
		// 	v.set_anchor(anchor);
		//
        //     if x > 0 && y > 0 {
        //         let vx = (x as f32 - 0.5) / max_dimension - 1.0;
        //         let vy = -1_f32*((y as f32 - 0.5) / max_dimension - 1.0);
		// 	    mesh.add_vertex([vx, vy, 0.0], [0_f32, 0_f32, 0_f32]);
        //     }
		//
		// }
		// println!("{}/{} vertices", mesh.vertices.len(), (width)*(height));
		// // Create triangles
		// for y in 0..height-1 {
		// 	for x in 0..width-1 {
        //         let i1 = (x + (width * y)) as Index;
        //         let i2 = (x + (width * y) + 1) as Index;
        //         let i3 = (x + (width * (y+1)) + 1) as Index;
        //         let i4 = (x + (width * (y+1))) as Index;
        //         if (x+y) & 1 == 1 {
        //             mesh.add_triangle(i1, i4, i2);
        //             mesh.add_triangle(i4, i3, i2);
        //         }
        //         else {
        //             mesh.add_triangle(i1, i3, i2);
        //             mesh.add_triangle(i4, i3, i1);
        //         }
		// 	}
		// }
		// println!("{}/{} triangles", mesh.triangles.len(), (width-1)*(height-1)*2);
		// println!("{}/{} edges", mesh.edges.len(), (width-1)*(height-1)*6);
		// mesh
	}
	pub fn from_image4(dynamic_image: DynamicImage) -> Self { 
        let offset = 0_f32;
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new();

		let max_dimension = (if width > height { width } else { height } as f32) / 2.0;
		// Create vertices
		for (x, y, pixel) in image.enumerate_pixels() {
			let v1 = mesh.add_vertex([(x as f32) / max_dimension - 1.0,        -1_f32*((y as f32) / max_dimension - 1.0), 0.0], pixel.0);
			let v2 = mesh.add_vertex([(x as f32 + 1.0-offset) / max_dimension - 1.0,  -1_f32*((y as f32) / max_dimension - 1.0), 0.0], pixel.0);
			let v3 = mesh.add_vertex([(x as f32) / max_dimension - 1.0,        -1_f32*((y as f32 + 1.0-offset) / max_dimension - 1.0), 0.0], pixel.0);
			let v4 = mesh.add_vertex([(x as f32 + 1.0-offset) / max_dimension - 1.0,  -1_f32*((y as f32 + 1.0-offset) / max_dimension - 1.0), 0.0], pixel.0);

            let xi = x as usize;
            let yi = y as usize;

            if xi == 0 {
                v1.set_anchor(ANCHOR_X);
                v3.set_anchor(ANCHOR_X);
            }
            if xi == width-1 {
                v2.set_anchor(ANCHOR_X);
                v4.set_anchor(ANCHOR_X);
            }
            if yi == 0 {
                v1.set_anchor(ANCHOR_Y);
                v2.set_anchor(ANCHOR_Y);
            }
            if yi == height-1 {
                v3.set_anchor(ANCHOR_Y);
                v4.set_anchor(ANCHOR_Y);
            }

            if xi == 0 && yi == 0 {
                v1.set_anchor(ANCHOR_BOTH);
            }
            if xi == width-1 && yi == 0 {
                v2.set_anchor(ANCHOR_BOTH);
            }
            if xi == 0 && yi == height-1 {
                v3.set_anchor(ANCHOR_BOTH);
            }
            if xi == width-1 && yi == height-1 {
                v4.set_anchor(ANCHOR_BOTH);
            } 
		}
		println!("{}/{} vertices", mesh.vertices.len(), 4*width*height);
		// Create triangles
        for x in 0..width {
            for y in 0..height {
                let i1 = 4*(x + (width * y)) as Index;

                mesh.add_triangle(i1, i1+2, i1+1);
                mesh.add_triangle(i1+2, i1+3, i1+1);

                if x < width-1 {
                    let i2 = 4*(x + (width * y) + 1) as Index;
                    // mesh.add_triangle(i1+1, i1+3, i2);
                    // mesh.add_triangle(i1+3, i2+2, i2);
                    mesh.add_triangle(i1+1, i2+2, i2);
                    mesh.add_triangle(i1+3, i2+2, i1+1);
                }
                if y < height-1 {
                    let i3 = 4*(x + (width * (y+1))) as Index;
                    // mesh.add_triangle(i1+2, i3, i1+3);
                    // mesh.add_triangle(i3, i3+1, i1+3);
                    mesh.add_triangle(i1+2, i3+1, i1+3);
                    mesh.add_triangle(i3, i3+1, i1+2);
                }
                if x < width-1 && y < height-1 {
                    let i2 = 4*(x + (width * y) + 1) as Index;
                    let i3 = 4*(x + (width * (y+1))) as Index;
                    let i4 = 4*(x + (width * (y+1)) + 1) as Index;
                    mesh.add_triangle(i1+3, i3+1, i2+2);
                    mesh.add_triangle(i3+1, i4, i2+2);
                }
			}
		}
        let tri_count = 2*width*height + 2*(width-1)*height + 2*width*(height-1) + 2*(width-1)*(height-1);
		println!("{}/{} triangles", mesh.triangles.len(), tri_count);
		println!("{}/{} edges", mesh.edges.len(), 3*tri_count);
		mesh
	}
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn triangle_count(&self) -> usize { self.triangles.len() }
	pub fn get_random_edge(&self) -> EdgePointer { 
		let i = rand::thread_rng().gen_range(0..self.edges.len()-1);
		self.edges[i].clone()
	}
	pub fn get_best_edge(&self) -> EdgePointer { 
        let mut index:usize = 0;
        let mut min: f32 = self.edges[0].energy();
        for (i,edge) in self.edges.iter().enumerate() {
            let energy = edge.energy();
            println!("{:?}", energy);
            if energy < min {
                min = energy;
                index = i;
            }
        }

		self.edges[index].clone()
	}
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
		let (v1, v2) = (edge.vertex(), edge.next().vertex());
		
		let (a1, a2) = (v1.anchor(), v2.anchor());
		let (p1, p2) = (v1.position(), v2.position());
		let (s1, s2) = (v1.state(), v2.state());
		let (c1, c2) = (v1.color(), v2.color());
		//println!("{:?}",v2.position());
        // println!("{:?},{} {:?},{}",p1,a1,p2,a2);
		let new_position = match (a1, a2) {
			(ANCHOR_NONE, ANCHOR_NONE) => average(p1, p2),
			(ANCHOR_BOTH, ANCHOR_NONE) => p1,
			(ANCHOR_NONE, ANCHOR_BOTH) => p2,

			(ANCHOR_BOTH, ANCHOR_X) if (p1[0] - p2[0]).abs() < EPSILON => p1,
			(ANCHOR_X, ANCHOR_BOTH) if (p1[0] - p2[0]).abs() < EPSILON => p2,
			(ANCHOR_X, ANCHOR_NONE) if (p1[0] - p2[0]).abs() < EPSILON  => p1,
			(ANCHOR_NONE, ANCHOR_X) if (p1[0] - p2[0]).abs() < EPSILON  => p2,

			(ANCHOR_BOTH, ANCHOR_Y) if (p1[1] - p2[1]).abs() < EPSILON => p1,
			(ANCHOR_Y, ANCHOR_BOTH) if (p1[1] - p2[1]).abs() < EPSILON => p2,
			(ANCHOR_Y, ANCHOR_NONE) if (p1[1] - p2[1]).abs() < EPSILON  => p1,
			(ANCHOR_NONE, ANCHOR_Y) if (p1[1] - p2[1]).abs() < EPSILON  => p2,

			(ANCHOR_X, ANCHOR_X) if (p1[0] - p2[0]).abs() < EPSILON  => average(p1, p2),
			(ANCHOR_Y, ANCHOR_Y) if (p1[1] - p2[1]).abs() < EPSILON  => average(p1, p2),
			_							=>	return Err("Invalid collapse".to_string()),
		};

		let new_color = average(c1, c2);
		
		let new_state = if s1 > s2 { s1 + 1 } else { s2 + 1 };
		// let new_anchor = match (a1, a2) {
		// q	(a, b) if a == b	=> 	a,
		// 	(a, ANCHOR_NONE)    =>	a,
		// 	(ANCHOR_NONE, b)    =>	b,
		// 	_ 					=>  panic!("How did you get here cotton eyed joe?"),
		// };
        let new_anchor = a1 | a2;
		
		// update the vertex
		v2.set_color(new_color);
		v2.set_position(new_position);
		v2.set_anchor(new_anchor);
		v2.set_state(new_state);
		//println!("> {:?}",v2.position());
		// update all the edges with the new vertex and collect all triangles that will be updated
		let mut triangles: HashSet<Index> = HashSet::new();
		let mut edges = Vec::new();
		
		for e in &mut edge.neighbors() {
			//e.set_vertex(&v_new);
			e.set_vertex(&v2);
			triangles.insert(e.triangle().index());
			edges.push(e.index());
		}
		for e in &mut edge.next().neighbors() {
			//e.set_vertex(&v_new);
			triangles.insert(e.triangle().index());
			//edges2.push(e.index());
		}
	
		// check if this results in the updated triangles overlaping any other triangles
		let t1 = edge.triangle().index();
		let t2 = if edge.has_opposite() { edge.opposite().unwrap().triangle().index() } else { t1 };
		assert!(triangles.remove(&t1));
		if edge.has_opposite() {
			assert!(triangles.remove(&t2)); // assertion fail
		}
		for modified_triangle_idx in &triangles {
			let mut modified_triangle = self.triangles[*modified_triangle_idx as usize].clone();
			modified_triangle.update();
			for triangle in &self.triangles {
				if triangles.contains(&triangle.index()) { continue; }
				if triangle.index() == t1 { continue; }
				if triangle.index() == t2 { continue; }
				if modified_triangle.overlaps(&triangle) {
					// OH NO REVERT THE CHANGES AND PRETEND THIS NEVER HAPPENED
					for e in edges {
						self.edges[e as usize].set_vertex(&v1);
					}
					v2.set_color(c2);
					v2.set_position(p2);
					v2.set_anchor(a2);
					v2.set_state(s2);
					//println!("< {:?}",v2.position());
					for modified_triangle_idx in &triangles {
						let modified_triangle = &mut self.triangles[*modified_triangle_idx as usize];
						modified_triangle.update();
					}
					return Err("Invalid collapse".to_string());
				}
			}
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
		
		// remove the triangles
		if edge.has_opposite() {
			assert_eq!(self.remove_triangle(edge.opposite().unwrap().triangle()), true);
		}
		assert_eq!(self.remove_triangle(edge.triangle()), true);
		
		// disable the vertex
		v1.set_state(0);
		// assert_eq!(self.remove_vertex(v1), true);
		// assert_eq!(self.remove_vertex(v2), true);
		
		// Return that we finished properly
		Ok(())
	}
    pub fn collapse_best_edge(&mut self) -> Result<(), String> {
        let mut edge_pq: BinaryHeap<_> = self.edges.iter().map(|x| x.clone()).map(Reverse).collect();
        while let Some(Reverse(edge)) = edge_pq.pop() {
            match self.collapse_edge(edge) {
                Ok(()) => return Ok(()),
                Err(e) => continue,
            }
        }
        Err("no valid collapses".to_string())
    }
    pub fn collapse_random_edge(&mut self) -> Result<(), String> {
        let mut indices: Vec<usize> = (0..self.edges.len()-1).collect();
        indices.shuffle(&mut rand::thread_rng());
        while let Some(i) = indices.pop() {
            let edge = self.edges[i].clone();
            match self.collapse_edge(edge) {
                Ok(()) => return Ok(()),
                Err(e) => continue,
            }
        }
        Err("no valid collapses".to_string())
    }
	
	fn add_vertex(&mut self, pos: Vec3f, color: Vec3f) -> VertexPointer { 
		self.vertices.push(VertexPointer::new(pos, color, self.vertices.len()));
		self.vertices.last().unwrap().clone()
	}
	fn add_edge(&mut self, start: &VertexPointer, end: &VertexPointer) -> EdgePointer {
        self.edges.push(EdgePointer::new(start, self.edges.len()));
        // self.energy_set.insert(self.edges.last().unwrap().clone());
		self.vertex_edge_map.insert((start.index(), end.index()), self.edges.last().unwrap().clone());
		self.edges.last().unwrap().clone()
	}
	fn add_triangle(&mut self, i1: Index, i2: Index, i3: Index) {
		let v1 = self.vertices[i1 as usize].clone();
		let v2 = self.vertices[i2 as usize].clone();
		let v3 = self.vertices[i3 as usize].clone();
		
		assert_eq!(self.find_edge(i1, i2), None);
		assert_eq!(self.find_edge(i2, i3), None);
		assert_eq!(self.find_edge(i3, i1), None);
		
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
		self.triangles.last_mut().unwrap().update();
		
		e1.set_triangle(self.triangles.last().unwrap());
		e2.set_triangle(self.triangles.last().unwrap());
		e3.set_triangle(self.triangles.last().unwrap());
	}
	
	fn remove_vertex(&mut self, vertex: VertexPointer) -> bool {
		let index = vertex.index() as usize;
		if index >= self.vertices.len() { return false; }
		
		self.vertices.swap_remove(index);
		
		if index == self.vertices.len() { return true; }
		assert_eq!(self.vertices[index].index() as usize, self.vertices.len());
		self.vertices[index].set_index(index as Index);
		
		true
	}
	fn remove_edge(&mut self, edge: EdgePointer) -> bool {
		let index = edge.index() as usize;
		if index > self.edges.len() { return false; }
		
		self.edges.swap_remove(index);
		if index == self.edges.len() { return true; }
		self.edges[index].set_index(index as Index);
		
		true
	}
	fn remove_triangle(&mut self, triangle: TrianglePointer) -> bool {
		let index = triangle.index() as usize;
		if index > self.triangles.len() { return false; }
		
		assert_eq!(self.remove_edge(triangle.edge().next().next()), true);
		assert_eq!(self.remove_edge(triangle.edge().next()), true);
		assert_eq!(self.remove_edge(triangle.edge()), true);
		
		self.triangles.swap_remove(index);
		if index == self.triangles.len() { return true; }
		self.triangles[index].set_index(index as Index);
		
		true
	}
}
