use std::cell::{Ref, RefCell};
use crate::vector::MyVector;
use rand::Rng;
use std::cmp::min;
use std::collections::HashMap;

const MIN_PARTICLE_SPEED: f32 = -2.0;
const MAX_PARTICLE_SPEED: f32 = 2.0;

const PRIME1: u64 = 6614058611;
const PRIME2: u64 = 7528850467;

#[derive(Debug)]
pub struct Particle {
    pub pos: MyVector,
    previous_pos: MyVector,
    velocity: MyVector,
    pub radius: f32,
    pub color: [u8; 4],
    id: u32,
}

impl Particle {
    pub fn get_display_pos(&self) -> (u32, u32) {
        (self.pos.x() as u32, self.pos.y() as u32)
    }

    fn update_position(&mut self, pos: MyVector) {
        self.previous_pos = self.pos;
        self.pos = pos;
    }

    fn get_border_x_min(&self) -> f32 {
        self.pos.x() - self.radius
    }

    fn get_border_x_max(&self) -> f32 {
        self.pos.x() + self.radius
    }
    fn get_border_y_min(&self) -> f32 {
        self.pos.y() - self.radius
    }

    fn get_border_y_max(&self) -> f32 {
        self.pos.y() + self.radius
    }

    fn set_color(&mut self, color: [u8; 4]) {
        self.color = color
    }
}
pub struct ParticleHashGrid {
    cell_size: u32,
    hash_grid_size: u64,
    hash_grid: HashMap<u64, Vec<u32>>,
}

impl ParticleHashGrid {
    pub fn new() -> Self {
        ParticleHashGrid {
            cell_size: 25,
            hash_grid_size: 10000000,
            hash_grid: HashMap::new(),
        }
    }

    fn cell_hash_from_pos(&self, particle_pos: &MyVector) -> u64 {
        let x = (particle_pos.x() as u32 / self.cell_size) as u64;
        let y = (particle_pos.y() as u32 / self.cell_size) as u64;
        let x_prime = x.checked_mul(PRIME1);
        let y_prime = y.checked_mul(PRIME2);
        if x_prime.is_some() && y_prime.is_some() {
            return (x_prime.unwrap() ^ y_prime.unwrap()) % self.hash_grid_size;
        }
        0u64
    }

    fn cell_hash_from_index(&self, x: i32, y: i32) -> u64 {
        let x_prime = (x as u64).checked_mul(PRIME1);
        let y_prime = (y as u64).checked_mul(PRIME2);
        if x_prime.is_some() && y_prime.is_some() {
            return (x_prime.unwrap() ^ y_prime.unwrap()) % self.hash_grid_size;
        }
        0u64
    }
    fn get_cell_particle_ids_from_hash(&self, hash: u64) -> Vec<u32> {
        let mut id_list = vec![];
        if let Some(l) = self.hash_grid.get(&hash) {
            l.iter().for_each(|id| {
                id_list.push(*id);
            });
        };
        id_list
    }

    fn clear_cell_grid(&mut self) {
        self.hash_grid.clear();
    }

    pub fn neighbour_search(&mut self, particles_map: &HashMap<u32, Particle>) {
        self.clear_cell_grid();
        self.push_cell_particle_ids(particles_map);
    }

    fn push_cell_particle_ids(&mut self, particles_map: &HashMap<u32, Particle>) {
        for (id, p) in particles_map.iter() {
            let hash = self.cell_hash_from_pos(&p.pos);
            if let Some(l) = self.hash_grid.get_mut(&hash) {
                l.push(*id);
            } else {
                let mut l = Vec::new();
                l.push(*id);
                self.hash_grid.insert(hash, l);
            }
        }
    }

    pub fn get_cell_particle_ids_from_pos(&self, pos: &MyVector) -> Vec<u32> {
        let hash = self.cell_hash_from_pos(pos);
        self.get_cell_particle_ids_from_hash(hash)
    }
}
pub struct ParticleWorld {
    velocity_damping: f32,
    collision_damping: f32,
    gravity: MyVector,
    rest_density: f32,
    k_near: f32,
    k: f32,
    interaction_radius: f32,
    width: u32,
    height: u32,
}

impl ParticleWorld {
    pub fn new(width: u32, height: u32) -> Self {
        ParticleWorld {
            velocity_damping: 1.0,
            collision_damping: 1.0,
            gravity: MyVector::new(0.0, 10.0),
            rest_density: 30.0,
            k_near: 3.0,
            k: 0.05,
            interaction_radius: 10.0,
            width,
            height,
        }
    }

    pub fn create_particles(&self, particles_map: &mut HashMap<u32, Particle>, amount: u32) {
        let x_particles = amount.isqrt();
        let y_particles = x_particles;
        let radius = 2;

        let mut x_start = self.width / 2 - x_particles * radius/ 2;
        let mut y_start = self.height / 2 - y_particles * radius / 2;
        let mut id = 0;
        for _ in 0..x_particles {
            for _ in 0..y_particles {
                let p = Particle {
                    pos: MyVector::new_u32(x_start, y_start),
                    previous_pos: MyVector::new_u32(x_start, y_start),
                    // velocity: self.get_random_speed() * 10.0,
                    velocity: MyVector::zero(),
                    radius: radius as f32,
                    color: [0xac, 0x00, 0xe6, 0xff],
                    id,
                };
                particles_map.insert(id, p);
                id += 1;
                x_start = min(x_start + radius, self.width);
            }
            x_start = self.width / 2 - x_particles * radius / 2;
            y_start = min(y_start + radius, self.height);
        }
    }
    pub fn predict_positions(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(id, p)| {
            let pos_delta = p.velocity * dt * self.velocity_damping;
            p.update_position(p.pos + pos_delta);
        });
    }

    pub fn compute_velocity(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(id, p)| {
            let velocity = (p.pos - p.previous_pos) * (1.0 / dt);
            p.velocity = velocity;
        })
    }

    pub fn apply_gravity(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(id, p)| {
            p.velocity = p.velocity + (self.gravity * dt);
        })
    }

    pub fn check_boundaries_gas(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(id, p)| {
            if p.get_border_x_min() <= 0.0 {
                p.pos.set_x(p.radius);
                p.velocity
                    .set_x(p.velocity.x() * -1.0 * self.collision_damping);
            } else if p.get_border_x_max() >= self.width as f32 {
                p.pos.set_x((self.width - 1) as f32 - p.radius);
                p.velocity
                    .set_x(p.velocity.x() * -1.0 * self.collision_damping);
            }
            if p.get_border_y_min() <= 0.0 {
                p.pos.set_y(p.radius);
                p.velocity
                    .set_y(p.velocity.y() * -1.0 * self.collision_damping);
            } else if p.get_border_y_max() >= self.height as f32 {
                p.pos.set_y((self.height - 1) as f32 - p.radius);
                p.velocity
                    .set_y(p.velocity.y() * -1.0 * self.collision_damping);
            }
        })
    }

    // no bounce
    pub fn check_boundaries_fluid(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(id, p)| {
            if p.get_border_x_min() <= 0.0 {
                p.pos.set_x(p.radius);
                p.update_position(p.pos);
            } else if p.get_border_x_max() >= self.width as f32 {
                p.pos.set_x((self.width - 1) as f32 - p.radius);
                p.update_position(p.pos);
            }
            if p.get_border_y_min() <= 0.0 {
                p.pos.set_y(p.radius);
                p.update_position(p.pos);
            } else if p.get_border_y_max() >= self.height as f32 {
                p.pos.set_y((self.height - 1) as f32 - p.radius);
                p.update_position(p.pos);
            }
        })
    }

    fn get_random_speed(&self) -> MyVector {
        MyVector::new(
            rand::rng().random_range(MIN_PARTICLE_SPEED..MAX_PARTICLE_SPEED),
            rand::rng().random_range(0.0..MAX_PARTICLE_SPEED),
        )
    }

    pub fn get_neighbour_cell_particle_ids_for_pos(
        &self,
        particles_map: &HashMap<u32, Particle>,
        hash_grid: &ParticleHashGrid,
        id: u32,
        id_list: &mut Vec<u32>,
        center_pos: &MyVector,
        influence_distance: f32,
    ) {
        let p = particles_map.get(&id).unwrap();
        let x = p.pos.x() as u32 / hash_grid.cell_size;
        let y = p.pos.y() as u32 / hash_grid.cell_size;

        let n_pos = (x as i32, y as i32);
        for pos in [
            (n_pos.0, n_pos.1),
            (n_pos.0, n_pos.1 + 1),
            (n_pos.0, n_pos.1 - 1),
            (n_pos.0 - 1, n_pos.1),
            (n_pos.0 + 1, n_pos.1),
            (n_pos.0 + 1, n_pos.1 + 1),
            (n_pos.0 - 1, n_pos.1 + 1),
            (n_pos.0 + 1, n_pos.1 - 1),
            (n_pos.0 - 1, n_pos.1 - 1),
        ]
        .iter()
        {
            let hash = hash_grid.cell_hash_from_index(pos.0, pos.1);
            for neighbour_id in hash_grid.get_cell_particle_ids_from_hash(hash) {
                let d = self.get_distance_of_particle_to_pos(particles_map, neighbour_id, center_pos);
                if d < influence_distance {
                    id_list.push(neighbour_id);
                }
            }
        }
    }

    fn get_neighbour_cell_particles<'a>(
        &self,
        hash_grid: &ParticleHashGrid,
        p: &Particle,
        particle_list: &mut Vec<u32>,
    ) {
        let x = p.pos.x() as u32 / hash_grid.cell_size;
        let y = p.pos.y() as u32 / hash_grid.cell_size;

        let n_pos = (x as i32, y as i32);
        for pos in [
            (n_pos.0, n_pos.1),
            (n_pos.0, n_pos.1 + 1),
            (n_pos.0, n_pos.1 - 1),
            (n_pos.0 - 1, n_pos.1),
            (n_pos.0 + 1, n_pos.1),
            (n_pos.0 + 1, n_pos.1 + 1),
            (n_pos.0 - 1, n_pos.1 + 1),
            (n_pos.0 + 1, n_pos.1 - 1),
            (n_pos.0 - 1, n_pos.1 - 1),
        ]
        .iter()
        {
            let hash = hash_grid.cell_hash_from_index(pos.0, pos.1);
            for neighbour_id in hash_grid.get_cell_particle_ids_from_hash(hash) {
                particle_list.push(neighbour_id);
            }
        }
    }

    fn get_distance_of_particle_to_pos(&self, particles_map: &HashMap<u32, Particle>, id: u32, pos: &MyVector) -> f32 {
        let p = particles_map.get(&id).unwrap();
        (p.pos - *pos).magnitude()
    }

    // fn get_density_of_particles(&self, p1: &Particle, p2: &Particle) -> (f32, f32) {
    //     let d = (p1.pos - p2.pos).magnitude();
    //     let q = d / self.interaction_radius;
    //     if q < 1.0 {
    //         return ((1.0 - q).powi(2), (1.0 - q).powi(3));
    //     }
    //     (0.0, 0.0)
    // }

    // fn get_displacement_of_particles(
    //     &self,
    //     dt: f32,
    //     d: MyVector,
    //     pressure: f32,
    //     pressure_near: f32,
    // ) -> MyVector {
    //     let q = d.magnitude() / self.interaction_radius;
    //     if q < 1.0 {
    //         d.normalize();
    //         let displacement =
    //             dt.powi(2) * (pressure * (1.0 - q) + pressure_near * (1.0 - q).powi(2));
    //         return d * displacement;
    //     }
    //     MyVector::zero()
    // }

    pub fn double_density_relaxiation(&self, particles_map: &mut HashMap<u32, Particle>, hash_grid: &ParticleHashGrid, dt: f32) {
        let id_list = particles_map.keys().map(|k| *k).collect::<Vec<_>>();

        for id in id_list.iter() {
            let mut density = 0.0;
            let mut density_near = 0.0;
            let mut neighbours = vec![];
            let p = particles_map.get(id).unwrap();
            let pos = p.pos;

            self.get_neighbour_cell_particles(hash_grid, &p, &mut neighbours);
            for neighbour_id in neighbours.iter() {
                if neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();

                let d = pos - neighbour.pos;
                let q = d.magnitude() / self.interaction_radius;
                if q < 1.0 {
                    density += (1.0 - q).powi(2);
                    density_near += (1.0 - q).powi(3);
                }
            }

            let pressure = self.k * (density - self.rest_density);
            let pressure_near = self.k_near * density_near;
            let mut this_displacement = MyVector::zero();

            for neighbour_id in neighbours.iter() {
                if neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();
                let mut d = neighbour.pos - pos;
                let q = d.magnitude() / self.interaction_radius;
                if q < 1.0 {
                    d = d.normalize();
                    let displacement_term =
                        dt.powi(2) * (pressure * (1.0 - q) + pressure_near * (1.0 - q).powi(2));
                    let displacement = d * displacement_term;
                    particles_map.entry(*neighbour_id).and_modify(|p| p.pos = p.pos + (displacement * 0.5));
                    this_displacement = this_displacement - (displacement * 0.5);
                }
            }
            particles_map.entry(*id).and_modify(|p| p.pos = p.pos + this_displacement);
        }
    }
}
