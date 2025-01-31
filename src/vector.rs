use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, SubAssign};

#[derive(Debug, Copy, Clone)]
pub struct MyVector {
    value: (f32, f32),
}

impl MyVector {
    pub fn new(x: f32, y: f32) -> Self {
        MyVector { value: (x, y) }
    }

    pub fn new_u32(x: u32, y: u32) -> Self {
        MyVector {
            value: (x as f32, y as f32),
        }
    }
    pub fn new_pair(xy: (f32, f32)) -> Self {
        MyVector { value: xy }
    }

    pub fn x(&self) -> f32 {
        self.value.0
    }

    pub fn set_x(&mut self, val: f32) {
        self.value.0 = val
    }

    pub fn y(&self) -> f32 {
        self.value.1
    }

    pub fn set_y(&mut self, val: f32) {
        self.value.1 = val
    }

    pub fn zero() -> Self {
        MyVector { value: (0.0, 0.0) }
    }

    pub fn is_zero(vector: &MyVector) -> bool {
        vector.value.0 == 0.0 && vector.value.1 == 0.0
    }

    pub fn magnitude(&self) -> f32 {
        self.magnitude_squared().sqrt()
    }

    pub fn magnitude_squared(&self) -> f32 {
        (self.value.0 * self.value.0 + self.value.1 * self.value.1)
    }
    pub fn dot(&self, rhs: Self) -> f32 {
        (self.value.0 * rhs.value.0) + (self.value.1 * rhs.value.1)
    }

    pub fn normalize(&self) -> Self {
        let magnitude = self.magnitude();
        if magnitude != 0.0 {
            return self.mul(1.0 / magnitude);
        }
        self.clone()
    }

    pub fn set_magnitude(&self, new_magnitude: f32) -> Self {
        self.normalize().mul(new_magnitude)
    }
}

impl Add for MyVector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        MyVector::new(self.value.0 + rhs.value.0, self.value.1 + rhs.value.1)
    }
}

impl AddAssign for MyVector {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for MyVector {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        MyVector::new(self.value.0 - rhs.value.0, self.value.1 - rhs.value.1)
    }
}

impl SubAssign for MyVector {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for MyVector {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        MyVector::new(self.value.0 * rhs.value.0, self.value.1 * rhs.value.1)
    }
}

impl MulAssign for MyVector {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Mul<f32> for MyVector {
    type Output = Self;

    fn mul(self, f: f32) -> Self::Output {
        MyVector::new(self.value.0 * f, self.value.1 * f)
    }
}
