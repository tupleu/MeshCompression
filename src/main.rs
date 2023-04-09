mod mesh;
use mesh::Mesh;
mod render;
use render::{Vertex, run};

const VERTICES: &[Vertex] = &[
    Vertex {
        position: [-0.5, -0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [0.0, 1.0, 0.0],
    },
    Vertex {
        position: [0.0, -0.5, 0.0],
        color: [1.0, 0.0, 0.0],
    },
    Vertex {
        position: [0.5, -0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [0.0, 0.0, 1.0],
    },
    Vertex {
        position: [0.25, 0.0, 0.0],
        color: [0.0, 1.0, 0.0],
    },
    Vertex {
        position: [0.0, 0.5, 0.0],
        // color: [0.0, 0.0, 0.0],
        color: [1.0, 0.0, 0.0],
    },
    Vertex {
        position: [-0.25, 0.0, 0.0],
        color: [0.0, 0.0, 1.0],
    },
    Vertex {
        position: [0.0, -0.225, 0.0],
        // color: [1.0, 1.0, 1.0],
        color: [0.0, 0.0, 0.0],
    },
];

const INDICES: &[u16] = &[
    0, 1, 5,
    1, 2, 3,
    3, 4, 5,
    1, 3, 6,
    3, 5, 6,
    5, 1, 6,
];

fn main() {
    pollster::block_on(run(VERTICES, INDICES));
}
