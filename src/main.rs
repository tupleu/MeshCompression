mod mesh;
use mesh::Mesh;

fn main() {
    println!("Hello, world!");

	let mut my_mesh = Mesh::new();
	let mut my_img_mesh = Mesh::from_image(8, 16);
	
	println!("{}", my_img_mesh.triangles.len());
	
}
