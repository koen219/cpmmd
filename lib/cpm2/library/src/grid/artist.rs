// For reading and opening files
use rayon::prelude::*;
use std::io::BufWriter;
use std::path::Path;
use std::{fs::File, sync::Arc};

use crate::graph::Graph;
use crate::grid::Grid2D;
use crate::Spin;

#[derive(Debug, Clone, Copy)]
pub struct Color(pub u8, pub u8, pub u8);

pub fn default_color_function(spin: Spin) -> Color {
    if spin == Spin(0) {
        return Color(0, 0, 0);
    }
    let colors = [
        Color(220, 20, 60),   // Crimson Red
        Color(0, 191, 255),   // Deep Sky Blue
        Color(50, 205, 50),   // Lime Green
        Color(255, 215, 0),   // Goldenrod Yellow
        Color(255, 105, 180), // Hot Pink
        Color(64, 224, 208),  // Turquoise
        Color(120, 81, 169),  // Royal Purple
        Color(255, 69, 0),    // Orange Red
        Color(173, 216, 230), // Light Blue
        Color(255, 140, 0),   // Dark Orange
    ];
    colors[((spin.0 - 1) as usize % colors.len()) as usize]
}

#[derive(Debug, Clone)]
pub struct Picture {
    data: Vec<u8>,
    sizex: usize,
    sizey: usize,
}

impl Picture {
    fn new(sizex: usize, sizey: usize) -> Picture {
        Picture {
            data: vec![0; sizex * sizey * 4],
            sizex,
            sizey,
        }
    }
    fn draw(&mut self, x: usize, y: usize, color: Color) {
        let index = x * 4 + y * 4 * self.sizex;
        self.data[index + 0] = color.0;
        self.data[index + 1] = color.1;
        self.data[index + 2] = color.2;
        self.data[index + 3] = 255;
    }
    pub fn scale(&mut self, scale: usize) {
        self.data = scale_image(&self.data, self.sizex, self.sizey, scale);
        self.sizex *= scale;
        self.sizey *= scale;
    }
}

fn scale_image(
    image: &[u8],
    original_width: usize,
    original_height: usize,
    scale: usize,
) -> Vec<u8> {
    let new_width = original_width * scale;
    let new_height = original_height * scale;

    let mut scaled_image = vec![0; new_width * new_height * 4]; // Assuming RGBA format

    for y in 0..new_height {
        for x in 0..new_width {
            // Find the corresponding pixel in the original image
            let orig_x = x / scale;
            let orig_y = y / scale;

            // Get the index in the original image
            let orig_idx = (orig_y * original_width + orig_x) * 4;

            // Get the index in the scaled image
            let new_idx = (y * new_width + x) * 4;

            // Copy RGB values
            scaled_image[new_idx] = image[orig_idx]; // Red
            scaled_image[new_idx + 1] = image[orig_idx + 1]; // Green
            scaled_image[new_idx + 2] = image[orig_idx + 2]; // Blue
            scaled_image[new_idx + 3] = image[orig_idx + 3]; // Blue
        }
    }

    scaled_image
}

// This function is slower than draw(). Maybe it can be sped up by doing less parallel threads?
pub fn parallel_draw<F>(grid: &Grid2D, color_function: F) -> Picture
where
    F: Fn(Spin) -> Color + Sync + Send + Copy,
{
    let sizex = grid.sizex;
    let sizey = grid.sizey;

    // Initialize the picture
    let mut picture = Picture::new(2 * sizex, 2 * sizey);
    let grid = Arc::new(grid);

    let points: Vec<(usize, usize)> = (0..sizex)
        .flat_map(|i| (0..sizey).map(move |j| (i, j)))
        .collect();

    let data: Vec<(usize, usize, Color)> = points
        .into_par_iter()
        .flat_map({
            move |(x, y)| {
                let (xx, yy) = (2 * x, 2 * y);
                let boundary_x = x + 1 < sizex;
                let boundary_y = y + 1 < sizey;
                let spin = grid.get(grid.from_positions(x, y));
                let color = color_function(spin);
                let black = Color(0, 0, 0);
                let mut pixels: Vec<(usize, usize, Color)> = Vec::with_capacity(4);
                pixels.push((xx, yy, color));
                if boundary_x && grid.get(grid.from_positions(x + 1, y)) == spin {
                    pixels.push((xx + 1, yy, color));
                } else {
                    pixels.push((xx + 1, yy, black));
                }
                if boundary_y && grid.get(grid.from_positions(x, y + 1)) == spin {
                    pixels.push((xx, yy + 1, color));
                } else {
                    pixels.push((xx, yy + 1, black));
                }
                if boundary_x && boundary_y && grid.get(grid.from_positions(x + 1, y + 1)) == spin {
                    pixels.push((xx + 1, yy + 1, color));
                } else {
                    pixels.push((xx + 1, yy + 1, black));
                }
                pixels
            }
        })
        .collect();
    data.into_iter()
        .for_each(|(x, y, color)| picture.draw(x, y, color));
    picture
}
pub fn draw<F>(grid: &Grid2D, color_function: F) -> Picture
where
    F: Fn(Spin) -> Color + Copy,
{
    let sizex = grid.sizex;
    let sizey = grid.sizey;

    // Initialize the picture
    let mut picture = Picture::new(2 * sizex, 2 * sizey);
    let grid = Arc::new(grid);

    (0..sizex)
        .flat_map(|i| (0..sizey).map(move |j| (i, j)))
        .into_iter()
        .flat_map({
            move |(x, y)| {
                let (xx, yy) = (2 * x, 2 * y);
                let boundary_x = x + 1 < sizex;
                let boundary_y = y + 1 < sizey;
                let spin = grid.get(grid.from_positions(x, y));
                let color = color_function(spin);
                let black = Color(0, 0, 0);
                let mut pixels: Vec<(usize, usize, Color)> = Vec::with_capacity(4);
                pixels.push((xx, yy, color));
                if boundary_x && grid.get(grid.from_positions(x + 1, y)) == spin {
                    pixels.push((xx + 1, yy, color));
                } else {
                    pixels.push((xx + 1, yy, black));
                }
                if boundary_y && grid.get(grid.from_positions(x, y + 1)) == spin {
                    pixels.push((xx, yy + 1, color));
                } else {
                    pixels.push((xx, yy + 1, black));
                }
                if boundary_x && boundary_y && grid.get(grid.from_positions(x + 1, y + 1)) == spin {
                    pixels.push((xx + 1, yy + 1, color));
                } else {
                    pixels.push((xx + 1, yy + 1, black));
                }
                pixels
            }
        })
        .for_each(|(x, y, color)| picture.draw(x, y, color));
    picture
}

pub fn draw_old<F>(grid: &Grid2D, color_function: F) -> Picture
where
    F: Fn(Spin) -> Color,
{
    let mut picture = Picture::new(2 * grid.sizex, 2 * grid.sizey);
    for x in 0..grid.sizex {
        for y in 0..grid.sizey {
            let spin = grid.get(grid.from_positions(x, y));
            let black = Color(0, 0, 0);
            let color = color_function(spin);
            let (xx, yy) = (2 * x, 2 * y);

            let x_in_bounds = x + 1 < grid.sizex;
            let y_in_bounds = y + 1 < grid.sizey;

            picture.draw(xx, yy, color);
            if x_in_bounds && grid.get(grid.from_positions(x + 1, y)) == spin {
                picture.draw(xx + 1, yy, color);
            } else {
                picture.draw(xx + 1, yy, black);
            }
            if y_in_bounds && grid.get(grid.from_positions(x, y + 1)) == spin {
                picture.draw(xx, yy + 1, color);
            } else {
                picture.draw(xx, yy + 1, black);
            }
            if x_in_bounds && y_in_bounds && grid.get(grid.from_positions(x + 1, y + 1)) == spin {
                picture.draw(xx + 1, yy + 1, color);
            } else {
                picture.draw(xx + 1, yy + 1, black);
            }
        }
    }
    picture
}

pub fn write_picture_as_png(picture: &Picture, path: &Path) {
    let width = picture.sizex as u32;
    let height = picture.sizey as u32;
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let mut encoder = png::Encoder::new(w, width, height); // Width is 2 pixels and height is 1.
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_compression(png::Compression::Fast);
    encoder.set_filter(png::FilterType::Sub);

    let mut writer = encoder.write_header().unwrap();
    writer
        .write_image_data(&picture.data)
        .expect("Output folder does not exist.") //.unwrap(); // Save
}
