pub mod artist;
pub mod ellipse;
mod grid2d;
mod grid3d;
mod initalize;

pub use artist as export2d;
pub use grid2d::Grid2D;
// pub use grid3d::export as export3d;
pub use grid3d::Grid3D;
pub use initalize::threedim as init3d;
pub use initalize::twodim as init2d;
