use bevy::prelude::Vec2;
use rand::Rng;
use std::collections::HashMap;
use std::vec::Vec;

const MIN_PARTICLE_SPEED: f32 = -2.0;
const MAX_PARTICLE_SPEED: f32 = 2.0;

const PRIME1: u64 = 6614058611;
const PRIME2: u64 = 7528850467;

#[derive(Debug)]
pub struct Particle {
    pub pos: Vec2,
    previous_pos: Vec2,
    pub velocity: Vec2,
    pub radius: f32,
    id: u32,
}

impl Particle {
    pub fn get_display_pos(&self) -> (u32, u32) {
        (self.pos.x as u32, self.pos.y as u32)
    }

    fn update_position(&mut self, pos: Vec2) {
        self.previous_pos = self.pos;
        self.pos = pos;
    }

    fn get_border_x_min(&self) -> f32 {
        self.pos.x - self.radius
    }

    fn get_border_x_max(&self) -> f32 {
        self.pos.x + self.radius
    }
    fn get_border_y_min(&self) -> f32 {
        self.pos.y - self.radius
    }

    fn get_border_y_max(&self) -> f32 {
        self.pos.y + self.radius
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
            cell_size: 20,
            hash_grid_size: 10000000,
            hash_grid: HashMap::new(),
        }
    }

    fn cell_hash_from_pos(&self, particle_pos: &Vec2) -> u64 {
        let x = (particle_pos.x as u32 / self.cell_size) as u64;
        let y = (particle_pos.y as u32 / self.cell_size) as u64;
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

    pub fn get_cell_particle_ids_from_pos(&self, pos: &Vec2) -> Vec<u32> {
        let hash = self.cell_hash_from_pos(pos);
        self.get_cell_particle_ids_from_hash(hash)
    }
}

#[derive(Debug, Copy, Clone)]
pub struct SimulationProperties {
    pub gravity: Vec2,
    pub rest_density: f32,
    pub k_near: f32,
    pub k: f32,
    pub interaction_radius: f32,
    pub sigma: f32,
    pub beta: f32,
    pub velocity_damping: f32,
}

impl SimulationProperties {
    pub fn new() -> Self {
        SimulationProperties {
            gravity: Vec2::new(0.0, 1.0),
            rest_density: 20.0,
            k_near: 4.0,
            k: 0.05,
            interaction_radius: 20.0,
            sigma: 0.0,
            beta: 0.02,
            velocity_damping: 1.0,
        }
    }
}

pub struct ParticleWorld {
    width: f32,
    height: f32,
    pub properties: SimulationProperties,
}

impl ParticleWorld {
    pub fn new(width: f32, height: f32) -> Self {
        ParticleWorld {
            width,
            height,
            properties: SimulationProperties::new(),
        }
    }

    pub fn change_simulation_properties(&mut self, properties: SimulationProperties) {
        self.properties = properties;
    }

    pub fn create_particles(
        &self,
        particles_map: &mut HashMap<u32, Particle>,
        amount: f32,
        radius: f32,
    ) {
        let x_particles = amount.sqrt();
        let y_particles = x_particles;

        let mut x_start = self.width / 2.0 - x_particles * radius / 2.0;
        let mut y_start = 0.0;
        let mut id = 0;
        let particle_spacing = radius;
        for _ in 0..x_particles as usize {
            for _ in 0..y_particles as usize {
                let p = Particle {
                    pos: Vec2::new(x_start, y_start),
                    previous_pos: Vec2::new(x_start, y_start),
                    // velocity: self.get_random_speed() * 10.0,
                    velocity: Vec2::ZERO,
                    radius,
                    id,
                };
                particles_map.insert(id, p);
                id += 1;
                x_start = f32::min(x_start + radius + particle_spacing, self.width);
            }
            x_start = self.width / 2.0 - x_particles * radius / 2.0;
            y_start = f32::min(y_start + radius + particle_spacing, self.height);
        }
    }

    pub fn reset_particles(&self, particles_map: &mut HashMap<u32, Particle>, radius: f32) {
        let x_particles = (particles_map.len() as f32).sqrt();
        let y_particles = x_particles;

        let mut x_start = self.width / 2.0 - x_particles * radius / 2.0;
        let mut y_start = 0.0;

        let particle_spacing = radius;
        let mut p_id = 0;
        for _ in 0..x_particles as usize {
            for _ in 0..y_particles as usize {
                particles_map.entry(p_id).and_modify(|p| {
                    p.pos = Vec2::new(x_start, y_start);
                    p.previous_pos = Vec2::new(x_start, y_start);
                    p.velocity = Vec2::ZERO;
                });

                p_id += 1;
                x_start = f32::min(x_start + radius + particle_spacing, self.width);
            }
            x_start = self.width / 2.0 - x_particles * radius / 2.0;
            y_start = f32::min(y_start + radius + particle_spacing, self.height);
        }
    }

    pub fn predict_positions(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            let pos_delta = p.velocity * dt * self.properties.velocity_damping;
            p.update_position(p.pos + pos_delta);
        });
    }

    pub fn compute_velocity(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            let velocity = (p.pos - p.previous_pos) * (1.0 / dt);
            p.velocity = velocity;
        })
    }

    pub fn apply_gravity(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            p.velocity = p.velocity + (self.properties.gravity * dt);
        })
    }

    // no bounce
    pub fn check_boundaries_fluid(&self, particles_map: &mut HashMap<u32, Particle>, _dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            if p.get_border_x_min() <= 0.0 {
                p.pos.x = p.radius;
                p.update_position(p.pos);
            } else if p.get_border_x_max() >= self.width as f32 {
                p.pos.x = (self.width - 1.0) - p.radius;
                p.update_position(p.pos);
            }
            if p.get_border_y_min() <= 0.0 {
                p.pos.y = p.radius;
                p.update_position(p.pos);
            } else if p.get_border_y_max() >= self.height as f32 {
                p.pos.y = (self.height - 1.0) - p.radius;
                p.update_position(p.pos);
            }
        })
    }

    fn get_random_speed(&self) -> Vec2 {
        Vec2::new(
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
        center_pos: &Vec2,
    ) {
        let p = particles_map.get(&id).unwrap();
        let x = p.pos.x as u32 / hash_grid.cell_size;
        let y = p.pos.y as u32 / hash_grid.cell_size;

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
                let d =
                    self.get_distance_of_particle_to_pos(particles_map, neighbour_id, center_pos);
                if d < self.properties.interaction_radius {
                    id_list.push(neighbour_id);
                }
            }
        }
    }

    fn get_neighbour_cell_particles(
        &self,
        hash_grid: &ParticleHashGrid,
        p: &Particle,
        particle_list: &mut Vec<u32>,
    ) {
        let x = p.pos.x as u32 / hash_grid.cell_size;
        let y = p.pos.y as u32 / hash_grid.cell_size;

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

    fn get_distance_of_particle_to_pos(
        &self,
        particles_map: &HashMap<u32, Particle>,
        id: u32,
        pos: &Vec2,
    ) -> f32 {
        let p = particles_map.get(&id).unwrap();
        (p.pos - *pos).length()
    }

    pub fn double_density_relaxiation(
        &self,
        particles_map: &mut HashMap<u32, Particle>,
        hash_grid: &ParticleHashGrid,
        dt: f32,
    ) {
        let dt_pow = dt.powi(2);

        for id in 0..particles_map.len() as u32 {
            let mut density = 0.0;
            let mut density_near = 0.0;
            let mut neighbours = vec![];
            let p = particles_map.get(&id).unwrap();
            let pos = p.pos;
            let mut neighbours_filtered: HashMap<u32, (Vec2, f32)> = HashMap::new();

            self.get_neighbour_cell_particles(hash_grid, &p, &mut neighbours);
            for neighbour_id in neighbours.iter() {
                if *neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();

                let d = neighbour.pos - pos;
                let q = d.length() / self.properties.interaction_radius;

                if q < 1.0 {
                    neighbours_filtered.insert(*neighbour_id, (d, q));
                    density += (1.0 - q).powi(2);
                    density_near += (1.0 - q).powi(3);
                }
            }

            let pressure = self.properties.k * (density - self.properties.rest_density);
            let pressure_near = self.properties.k_near * density_near;
            let mut this_displacement = Vec2::ZERO;

            for (neighbour_id, (d, q)) in neighbours_filtered.iter() {
                let d = d.normalize();
                let displacement_term =
                    dt_pow * (pressure * (1.0 - q) + pressure_near * (1.0 - q).powi(2));
                let displacement = d * displacement_term;
                particles_map
                    .entry(*neighbour_id)
                    .and_modify(|p| p.pos = p.pos + (displacement * 0.5));
                this_displacement = this_displacement - (displacement * 0.5);
            }
            particles_map
                .entry(id)
                .and_modify(|p| p.pos = p.pos + this_displacement);
        }
    }

    pub fn viscosity(
        &self,
        particles_map: &mut HashMap<u32, Particle>,
        hash_grid: &ParticleHashGrid,
        dt: f32,
    ) {
        for id in 0..particles_map.len() as u32 {
            let mut neighbours = vec![];
            let p = particles_map.get(&id).unwrap();
            let p_pos = p.pos;
            let p_velocity = p.velocity;

            self.get_neighbour_cell_particles(hash_grid, &p, &mut neighbours);
            for neighbour_id in neighbours.iter() {
                if *neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();
                let d = neighbour.pos - p_pos;
                let q = d.length() / self.properties.interaction_radius;
                let v_diff = p_velocity - neighbour.velocity;

                if q < 1.0 {
                    let d = d.normalize();
                    let u = v_diff.dot(d);

                    if u > 0.0 {
                        let i_term = dt
                            * (1.0 - q)
                            * (self.properties.sigma * u + self.properties.beta * u * u);
                        let i = d * i_term;

                        particles_map
                            .entry(id)
                            .and_modify(|p| p.velocity = p.velocity - (i * 0.5));
                        particles_map
                            .entry(*neighbour_id)
                            .and_modify(|p| p.velocity = p.velocity + (i * 0.5));
                    }
                }
            }
        }
    }
}
