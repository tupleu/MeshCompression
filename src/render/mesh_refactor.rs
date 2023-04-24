use std::rc::{Rc, Weak};
use std::cell::RefCell;
use image::{DynamicImage, Pixel};
use std::collections::{HashMap, HashSet};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::cmp;
use kdtree::{KdTree, distance};

use ordered_float::OrderedFloat;
use indicatif::{ProgressBar, ProgressStyle, ProgressState, HumanDuration};
use std::{fmt::Write};
use std::time::Duration;


type Index = u32;
type Vec3f = [f32; 3];
type Vec2 = [u32; 2];

const ANCHOR_NONE: u16 = 0;
const ANCHOR_X: u16 = 1;
const ANCHOR_Y: u16 = 2;
const ANCHOR_BOTH: u16 = 3;

pub const RANDOM: usize = 0;

struct Square {
	pub upper_left: Vec2,
	pub upper_right: Vec2,
	pub lower_left: Vec2,
	pub lower_right: Vec2,
	pub pixel: Vec3f,
}

impl Square {
	pub fn contains(&self, point: Vec2) -> bool {
		point[0] >= self.x() && 
			point[0] <= self.x() + self.width() &&
			point[1] >= self.y() &&
			point[1] <= self.y() + self.height()
	}
	pub fn height(&self) -> u32 {
		self.upper_right[1] - self.lower_right[1]
	}
	pub fn width(&self) -> u32 {
		self.upper_right[0] - self.upper_left[0]
	}
	pub fn x(&self) -> u32 {
		self.lower_left[0]
	}
	pub fn y(&self) -> u32 {
		self.lower_left[1]
	}
	pub fn points(&self) -> [Vec2; 4] {
		[self.lower_left, self.lower_right, self.upper_right, self.upper_left]
	}
}

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

fn new_progress_bar(length: u64) -> ProgressBar {
	let bar = ProgressBar::new(length);
	bar.set_style(ProgressStyle::with_template("[{elapsed_precise}] [{wide_bar}] {percent}% ({smoothed_eta})")
		.unwrap()
		.with_key(
		  "smoothed_eta",
		  |s: &ProgressState, w: &mut dyn Write| match (s.pos(), s.len()) {
			  (pos, Some(len)) => write!(
				  w,
				  "{:#}",
				  HumanDuration(Duration::from_millis(
					  (s.elapsed().as_millis() * (len as u128 - pos as u128) / (pos as u128))
						  as u64
				  ))
			  )
			  .unwrap(),
			  _ => write!(w, "-").unwrap(),
		  })
		.progress_chars("██.."));
	bar
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

#[derive(Debug, Clone, PartialEq)]
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
	pub fn neighboors(&self) -> Vec<EdgePointer> {
		let self_vi = self.vertex().index();
		let mut results = Vec::new();
		let mut current_edge = self.clone();
		let mut previous_edge_idx = current_edge.index();
		loop {
			results.push(current_edge.clone());
			previous_edge_idx = current_edge.index();
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap().next();
			if current_edge.index() == self.index() { break; }
			assert_eq!(current_edge.vertex().index(), self_vi);
		}
		current_edge = self.next().next();
		if !current_edge.has_opposite() { return results; }
		current_edge = current_edge.opposite().unwrap();
		loop {
			if previous_edge_idx == current_edge.index() { break; }
			results.insert(0, current_edge.clone());
			current_edge = current_edge.next().next();
			if !current_edge.has_opposite() { break; }
			current_edge = current_edge.opposite().unwrap();
			if current_edge.index() == self.index() { break; }
			assert_eq!(current_edge.vertex().index(), self_vi);
		}
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
	pub fn update(&mut self, map: &mut KdTree<f32, Index, Vec3f>) {
		eprintln!("{:?}", &self.centroid());
		map.add(self.centroid(), self.index());
		match map.remove(&self.centroid(), &self.index()) {
			Err(e) => panic!("{}", e),
			Ok(_) => {},
		};
		self.update_centroid();
		map.add(self.centroid(), self.index());
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
	}
	pub fn overlaps(&self, triangle: &TrianglePointer) -> bool {
		// do a quick check for circular overlap
		let d = distance(self.centroid(), triangle.centroid());
		if d > self.circumscribed_radius() + triangle.circumscribed_radius() {
			return false;
		}
		// Check each of the edges
		for edge1 in self.edges() {
            for edge2 in triangle.edges() {
                if edge1.intersects(&edge2) {
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
	vertex_edge_map: HashMap<(Index, Index), EdgePointer>,
	history: Vec<(Index, Option<Index>, Index, Index, (f32, f32))>, // vl, vr, vs, (dx dy), (dx, dy)
	r: StdRng,
	triangle_map: KdTree<f32, Index, Vec3f>,
	history_map: KdTree<f32, Index, Vec3f>,
}

impl Mesh {
	// Constructors ///////////////////////////////////////////////////////////////////////////
	fn new(r: u64) -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            triangles: Vec::new(),
			vertex_edge_map: HashMap::new(),
			history: Vec::new(), 
			triangle_map: KdTree::new(3),
			history_map: KdTree::new(3),
			//KdTree2::build(vec![]),
			r: StdRng::seed_from_u64(r),
        }
    }
	
	pub fn from_image(dynamic_image: DynamicImage, r: u64) -> Self {
		let width = dynamic_image.width() as usize;
		let height = dynamic_image.height() as usize;
		println!("{}x{}", width, height);
		let image = dynamic_image.to_rgb32f();
		
		let mut mesh = Mesh::new(r);
		
		
		let max_dimension = (if width > height { width - 1 } else { height - 1 } as f32) / 2.0;

		// Create vertices
		println!("Creating verticies...");
		let mut bar = new_progress_bar((width * height) as u64);
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
			bar.inc(1);
		}
		bar.finish();
		println!("Made {:?} verticies!", mesh.vertices.len());
		
		// Create triangles
		println!("Creating triangles and edges...");
		bar = new_progress_bar(((width-1)*(height-1)) as u64);

		for y in 0..height-1 {
			for x in 0..width-1 {
				let i1 = (x + (width * y)) as Index;
				let i2 = (x + (width * y) + 1) as Index;
				let i3 = (x + (width * (y+1)) + 1) as Index;
				let i4 = (x + (width * (y+1))) as Index;
				let t1 = mesh.add_triangle(i1, i4, i2, true);
				let t2 = mesh.add_triangle(i4, i3, i2, true);
				bar.inc(1);
			}
		}
		bar.finish();
		println!("Made {:?} edges!", mesh.edges.len());
		println!("Made {:?} triangles!", mesh.triangles.len());
		//assert_eq!(mesh.edges.len(), (width-1)*(height-1)*6);
		mesh
	}
	
	// Getters ////////////////////////////////////////////////////////////////////////////////
	pub fn triangle_count(&self) -> usize { self.triangles.len() }
	pub fn edge_count(&self) -> usize { self.edges.len() }
	pub fn get_random_edge(&mut self) -> EdgePointer {
		let i = self.r.gen_range(0..self.edges.len()-1);
		self.edges[i].clone()
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
	pub fn undo_n_edge_collapses_at(&mut self, n: usize, at: Vec3f) {
		println!("Uncollapsing {} edges...", n);
		if n == 0 { return; }
		let bar = new_progress_bar(n as u64);
		bar.inc(1);
		for _ in 0..n {
			if self.history.len() == 0 { break; }
			let mut history_index = 0;
			let mut found = false;
			for (point, index) in self.history_map.iter_nearest(&at, &distance::squared_euclidean).unwrap() {
				bar.tick();
				let (il, ir, is, _, (_, _)) = self.history[*index as usize];
				let (vl, vs) = ( 
								self.vertices[il as usize].clone(),
								self.vertices[is as usize].clone(),
						   );
				let vr = if ir.is_some() { Some(self.vertices[ir.unwrap() as usize].clone()) } else { None };
				// keep looking if the edge cannot be uncollapsed
				if vl.state() == 0 || (vr.is_some() && vr.unwrap().state() == 0) || vs.state() == 0 { println!("oof"); continue; }
				history_index = *index;
				found = true;
				break;
			};
			if !found { continue; }
			self.undo_edge_collapse(history_index as usize);
			bar.inc(1);
		}
		bar.finish();
	}
	
	pub fn do_n_edge_collapses(&mut self, n: usize, method: usize) {
		let mut m = n;
		println!("Collapsing {} edges...", n);
		if n == 0 { return; }
		let bar = new_progress_bar(n as u64);
		bar.inc(1);
		while m > 0 {
			bar.tick();
			println!("{}", self.edges.len());
			let edge = match method {
				_ => self.get_random_edge(),
			};
			match self.collapse_edge(edge) {
                Err(e) => continue,
				_ => {},
            }
			
			bar.inc(1);
			m -= 1;
		}
		bar.finish();
	}
	
	pub fn undo_last_edge_collapse(&mut self) {
		if self.history.len() == 0 { return; }
		self.undo_edge_collapse(self.history.len() - 1);
	}
	fn undo_edge_collapse(&mut self, history_idx: usize) -> bool {
		let (il, ir, is, it, (dx, dy)) = self.history[history_idx];
		// grap the verticies
		let (vl, vs, vt) = ( 
								self.vertices[il as usize].clone(),
								self.vertices[is as usize].clone(),
								self.vertices[it as usize].clone()
						   );
		let vr = if ir.is_some() { Some(self.vertices[ir.unwrap() as usize].clone()) } else { None };
		let has_r = vr.is_some();
		// stop if the edge cannot be uncollapsed
		if vl.state() == 0 || (has_r && vr.unwrap().state() == 0) || vs.state() == 0 { return false; }
		// get any associated edge
		let mut edge_check = None;
		for e in &self.edges {
			if e.vertex().index() == is {
				edge_check = Some(e);
				break;
			}
		}
		if edge_check.is_none() { return false; }
		let edge = edge_check.unwrap();
		// calculate the old position
		let p_cur = vs.pos();
		let p_old = [vs.pos()[0] + dx, vs.pos()[1] + dy, 0.0];
		// get the neighborhood
		let mut edges = edge.neighboors();
		// find what edges left and right are
		let mut left_edge_idx = 0;
		let mut right_edge_idx = 0;
		let mut i = 0;
		for edge in &edges {
			i += 1;
			if !edge.has_opposite() { continue; }
			if edge.opposite().unwrap().vertex().index() == il {
				left_edge_idx = i - 1;
			}
			else if Some(edge.opposite().unwrap().vertex().index()) == ir {
				right_edge_idx = i - 1;
			}
		}
		// reorder these so that the left edge is at index 0
		while left_edge_idx != 0 {
			let temp_i = edges.remove(0);
			edges.push(temp_i);
			left_edge_idx -= 1;
			if right_edge_idx == 0 {
				right_edge_idx = edges.len();
			}
			right_edge_idx -= 1;
		}
		//assert_eq!("{}", edges[left_edge_idx].opposite().unwrap().vertex().index(), ir);
		// get the ordering
		let indicies = [il, ir.unwrap_or(0), is, it];
		let mut order = [1, 1, 1, 1];
		for n in 0..4 { for m in 0..4 {
			if n == m { continue; }
			if indicies[n] > indicies[m] { order[n as usize] = order[n as usize] + 1; }
		}}
		//println!("{:?}", order);
		// move the edges to their proper positions
		let len = edges.len();
		for n in 0..edges.len() {
			if n < right_edge_idx { continue; }
			edges[n].set_vertex(&vt);
		}
		if has_r { edges[right_edge_idx].set_vertex(&vs); }
		edges[left_edge_idx].set_vertex(&vt);
		// oh yea, move this vertex as well
		vs.set_position(p_old);
		vt.set_state(vs.state());
		// re-add the triangle(s) and connect opposites
		let tl = self.add_triangle(il, it, is, false);
		if has_r {
			let tr = self.add_triangle(ir.unwrap(), is, it, false);
			tr.edge().next().set_opposite(&Some(tl.edge().next()));
			tl.edge().next().set_opposite(&Some(tr.edge().next()));
			
			edges[right_edge_idx].opposite().unwrap().set_opposite(&Some(tr.edge().next().next()));
			tr.edge().next().next().set_opposite(&edges[right_edge_idx].opposite());
			
			edges[right_edge_idx].set_opposite(&Some(tr.edge()));
			tr.edge().set_opposite(&Some(edges[right_edge_idx].clone()));
		}
		edges[left_edge_idx].opposite().unwrap().set_opposite(&Some(tl.edge().next().next()));
		tl.edge().next().next().set_opposite(&edges[left_edge_idx].opposite());
		
		edges[left_edge_idx].set_opposite(&Some(tl.edge()));
		tl.edge().set_opposite(&Some(edges[left_edge_idx].clone()));
		// done!
		self.history.pop();
		
		true
	}
	pub fn collapse_edge(&mut self, edge: EdgePointer) -> Result<(), String> {
		let (v1, v2) = (edge.vertex(), edge.next().vertex());
		
		let (a1, a2) = (v1.anchor(), v2.anchor());
		let (p1, p2) = (v1.position(), v2.position());
		let (s1, s2) = (v1.state(), v2.state());
		let (c1, c2) = (v1.color(), v2.color());

		let safe = s1 == s2; 

		let new_position = match (a1, a2) {
			(ANCHOR_NONE, ANCHOR_NONE) 	=>	average(p1, p2),
			(ANCHOR_BOTH, ANCHOR_NONE) 	=>	p1,
			(ANCHOR_NONE, ANCHOR_BOTH) 	=>	p2,
			_							=>	return Err("Invalid collapse".to_string()),
		};
		let new_color = average(c1, c2);
		let new_state = if s1 > s2 { s1 + 1 } else { s2 + 1 };
		let new_anchor = match (a1, a2) {
			(a, b) if a == b			=> 	a,
			(a, b) if b == ANCHOR_NONE 	=>	a,
			(a, b) if a == ANCHOR_NONE 	=>	b,
			_ 							=>  panic!("How did you get here cotton eyed joe?"),
		};
		// update the vertex
		v2.set_color(new_color);
		v2.set_position(new_position);
		v2.set_anchor(new_anchor);
		v2.set_state(new_state);
		// update all the edges with the new vertex and collect all triangles that will be updated
		let mut triangles: HashSet<Index> = HashSet::new();
		let mut edges = Vec::new();
		for e in &mut edge.neighboors() {
			e.set_vertex(&v2);
			triangles.insert(e.triangle().index());
			edges.push(e.index());
		}
		if edges.len() == 0 { return Err("Invalid collapse".to_string()); }
		
		if false {
			for e in &mut edge.next().neighboors() {
				triangles.insert(e.triangle().index());
			}
			// check if this results in the updated triangles overlaping any other triangles
			let t1 = edge.triangle().index();
			let t2 = if edge.has_opposite() { edge.opposite().unwrap().triangle().index() } else { t1 };
			assert!(triangles.remove(&t1));
			if edge.has_opposite() {
				assert!(triangles.remove(&t2)); // assertion fail
			}
			// update the triangles and get the center centroid
			let mut centroid = [0.0,0.0,0.0];
			for modified_triangle_idx in &triangles {
				let mut modified_triangle = self.triangles[*modified_triangle_idx as usize].clone();
				modified_triangle.update(&mut self.triangle_map);
				centroid[0] += modified_triangle.centroid()[0];
				centroid[1] += modified_triangle.centroid()[1];
				centroid[2] += modified_triangle.centroid()[2];
			}
			centroid[0] /= triangles.len() as f32;
			centroid[1] /= triangles.len() as f32;
			centroid[2] /= triangles.len() as f32;
			// get the circumscribed radius that encircles all triangles and the triangles they may overlap
			let mut radius = 0.0;
			for modified_triangle_idx in &triangles {
				let triangle = &self.triangles[*modified_triangle_idx as usize];
				let radius_new = distance(triangle.centroid(), [0.0, 0.0, 0.0]) + triangle.circumscribed_radius();
				if radius_new > radius { radius = radius_new; }
			}
			
			// gather the possible triangles that may be overlapable
			let overlaping_triangles_idx = self.triangle_map.within(&centroid, radius, &distance::squared_euclidean).unwrap();
			// make sure no overlap occurs
			for modified_triangle_idx in &triangles {
				
				let modified_triangle = self.triangles[*modified_triangle_idx as usize].clone();
				for triangle_idx in &overlaping_triangles_idx {
					//println!("?");
					if triangles.contains(triangle_idx.1) { continue; }
					if *triangle_idx.1 == t1 { continue; }
					if *triangle_idx.1 == t2 { continue; }
					if modified_triangle.overlaps(&self.triangles[*triangle_idx.1 as usize]) {
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
							modified_triangle.update(&mut self.triangle_map);
						}
						return Err("Invalid collapse".to_string());
					}
				}
			}
		}
		// update the history
		
		self.history.push(( edge.next().next().vertex().index(),
							if edge.has_opposite() { Some(edge.opposite().unwrap().next().next().vertex().index()) } else { None },
							v2.index(),
							v1.index(),
							(p2[0] - new_position[0], p2[1] - new_position[1])
						));
		self.history_map.add(v1.pos(), self.history.len() as u32);
		//old = dx + new
		//old - new = dx
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
				op_edge.next().opposite().unwrap().set_opposite(&op_edge.next().next().opposite());
			}
			if op_edge.next().next().has_opposite() {
				op_edge.next().next().opposite().unwrap().set_opposite(&op_edge.next().opposite());
			}
		}
		// remove the triangles
		if edge.has_opposite() {
			assert_eq!(self.remove_triangle(edge.opposite().unwrap().triangle()), true);
		}
		assert_eq!(self.remove_triangle(edge.triangle()), true);
		
		// disable the vertex
		v1.set_state(0);
		
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
	fn add_triangle(&mut self, i1: Index, i2: Index, i3: Index, connect: bool) -> TrianglePointer {
		let v1 = self.vertices[i1 as usize].clone();
		let v2 = self.vertices[i2 as usize].clone();
		let v3 = self.vertices[i3 as usize].clone();
		if connect {
			//assert_eq!(self.find_edge(i1, i2), None);
			//assert_eq!(self.find_edge(i2, i3), None);
			//assert_eq!(self.find_edge(i3, i1), None);
		}
		let mut e1 = self.add_edge(&v1, &v2);
		let mut e2 = self.add_edge(&v2, &v3);
		let mut e3 = self.add_edge(&v3, &v1);
		
		e1.set_next(&e2);
		e2.set_next(&e3);
		e3.set_next(&e1);
		
		if connect {
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
		}

		self.triangles.push(TrianglePointer::new(&e1, self.triangles.len()));
		let triangle = self.triangles.last_mut().unwrap();
		triangle.update_centroid();
		triangle.update_circumscribed_radius();
		self.triangle_map.add(triangle.centroid(), triangle.index());
		
		e1.set_triangle(triangle);
		e2.set_triangle(triangle);
		e3.set_triangle(triangle);
		
		triangle.clone()
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
		
		self.triangle_map.remove(&triangle.centroid(), &triangle.index());
		self.triangles.swap_remove(index);
		
		if index == self.triangles.len() { return true; }
		
		let moved_triangle = &mut self.triangles[index];
		self.triangle_map.remove(&moved_triangle.centroid(), &moved_triangle.index());
		moved_triangle.set_index(index as Index);
		self.triangle_map.add(moved_triangle.centroid(), moved_triangle.index());
		
		
		true
	}
}
