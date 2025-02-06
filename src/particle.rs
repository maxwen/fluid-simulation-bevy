use bevy::prelude::Vec2;
use rand::Rng;
use std::collections::HashMap;
use std::vec::Vec;
use serde::{Deserialize, Serialize};

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

#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct ParticleWorldProperties {
    pub gravity: Vec2,
    pub rest_density: f32,
    pub k_near: f32,
    pub k: f32,
    pub interaction_radius: f32,
    pub sigma: f32,
    pub beta: f32,
    pub particle_num: f32,
    pub particle_radius: f32,
}

impl ParticleWorldProperties {
    pub fn new() -> Self {
        ParticleWorldProperties {
            gravity: Vec2::new(0.0, 1.0),
            rest_density: 30.0,
            k_near: 4.0,
            k: 0.03,
            interaction_radius: 20.0,
            sigma: 0.0,
            beta: 0.03,
            particle_num: 2000.0,
            particle_radius: 4.0,
        }
    }
}

pub struct CircleShape {
    pos: Vec2,
    radius: f32,
    active: bool,
    dampening_factor: f32,
    include: bool,
}

impl CircleShape {
    pub fn new(pos: Vec2, radius: f32, active: bool) -> Self {
        CircleShape {
            pos,
            radius,
            active,
            dampening_factor: 0.9,
            include: false,
        }
    }

    fn get_direction_out(&self, pos: Vec2) -> Vec2 {
        if self.active {
            let radius_squared = self.radius * self.radius;
            let d = pos - self.pos;
            if self.include {
                let d_out = radius_squared + self.radius / 2.0;
                let d_in = radius_squared - self.radius / 2.0;
                if d.length_squared() < d_in {
                    let p = self.dampening_factor;
                    let d = d.normalize();
                    return d * -p;
                } else if d.length_squared() < d_out {
                    let p = self.dampening_factor;
                    let d = d.normalize();
                    return d * p;
                }
            } else {
                if d.length_squared() < radius_squared {
                    let p = self.dampening_factor * 2.0;
                    let d = d.normalize();
                    return d * p;
                }
            }
        }
        Vec2::ZERO
    }

    pub fn set_position(&mut self, pos: Vec2) {
        self.pos = pos;
        self.active = true;
    }

    pub fn disable(&mut self) {
        self.active = false;
    }

    pub fn enable(&mut self, include: bool) {
        self.active = true;
        self.include = include;
    }

    pub fn is_active(&self) -> bool {
        self.active
    }
}

pub struct ParticleWorld {
    width: f32,
    height: f32,
    properties: ParticleWorldProperties,
    grid: ParticleHashGrid
}

impl ParticleWorld {
    pub fn new(width: f32, height: f32, properties: ParticleWorldProperties) -> Self {
        ParticleWorld {
            width,
            height,
            properties,
            grid: ParticleHashGrid::new()
        }
    }

    pub fn change_simulation_properties(&mut self, properties: ParticleWorldProperties) {
        self.properties = properties;
    }

    pub fn create_particles(&self, particles_map: &mut HashMap<u32, Particle>) {
        particles_map.clear();
        let x_particles = self.properties.particle_num.sqrt();
        let y_particles = x_particles;

        let mut x_start = self.width / 2.0 - x_particles * self.properties.particle_radius / 2.0;
        let mut y_start = 0.0;
        let mut id = 0;
        let particle_spacing = self.properties.particle_radius;
        for _ in 0..x_particles as usize {
            for _ in 0..y_particles as usize {
                let p = Particle {
                    pos: Vec2::new(x_start, y_start),
                    previous_pos: Vec2::new(x_start, y_start),
                    // velocity: self.get_random_speed() * 10.0,
                    velocity: Vec2::ZERO,
                    radius: self.properties.particle_radius,
                    id,
                };
                particles_map.insert(id, p);
                id += 1;
                x_start = f32::min(
                    x_start + self.properties.particle_radius + particle_spacing,
                    self.width,
                );
            }
            x_start = self.width / 2.0 - x_particles * self.properties.particle_radius / 2.0;
            y_start = f32::min(
                y_start + self.properties.particle_radius + particle_spacing,
                self.height,
            );
        }
    }

    pub fn predict_positions(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            let pos_delta = p.velocity * dt;
            p.update_position(p.pos + pos_delta);
        });
    }

    pub fn compute_velocity(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            let velocity = (p.pos - p.previous_pos) * (1.0 / dt);
            p.velocity = velocity;
        })
    }

    pub fn compute_velocity_and_check_bounds(&self, particles_map: &mut HashMap<u32, Particle>, dt: f32) {
        particles_map.iter_mut().for_each(|(_id, p)| {
            if p.get_border_x_min() <= 0.0 {
                p.pos.x = p.radius;
                p.update_position(p.pos);
            } else if p.get_border_x_max() >= self.width {
                p.pos.x = (self.width - 1.0) - p.radius;
                p.update_position(p.pos);
            }
            if p.get_border_y_min() <= 0.0 {
                p.pos.y = p.radius;
                p.update_position(p.pos);
            } else if p.get_border_y_max() >= self.height {
                p.pos.y = (self.height - 1.0) - p.radius;
                p.update_position(p.pos);
            }

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
            } else if p.get_border_x_max() >= self.width {
                p.pos.x = (self.width - 1.0) - p.radius;
                p.update_position(p.pos);
            }
            if p.get_border_y_min() <= 0.0 {
                p.pos.y = p.radius;
                p.update_position(p.pos);
            } else if p.get_border_y_max() >= self.height {
                p.pos.y = (self.height - 1.0) - p.radius;
                p.update_position(p.pos);
            }
        })
    }

    fn check_shape_collision(
        &self,
        particles_map: &mut HashMap<u32, Particle>,
        circle_shape: &CircleShape,
        _dt: f32,
    ) {
        if circle_shape.active {
            particles_map.iter_mut().for_each(|(_id, p)| {
                let d = circle_shape.get_direction_out(p.pos);
                if d != Vec2::ZERO {
                    p.pos = p.pos + d;
                }
            });
        }
    }

    #[allow(dead_code)]
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
        let interaction_radius_squared =
            self.properties.interaction_radius * self.properties.interaction_radius;

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
                let d = self.get_distance_squared_of_particle_to_pos(
                    particles_map,
                    neighbour_id,
                    center_pos,
                );
                if d < interaction_radius_squared {
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

    fn get_distance_squared_of_particle_to_pos(
        &self,
        particles_map: &HashMap<u32, Particle>,
        id: u32,
        pos: &Vec2,
    ) -> f32 {
        let p = particles_map.get(&id).unwrap();
        (p.pos - *pos).length_squared()
    }

    pub fn double_density_relaxiation(
        &self,
        particles_map: &mut HashMap<u32, Particle>,
        dt: f32,
    ) {
        let dt_pow = dt.powi(2);
        let interaction_radius_squared =
            self.properties.interaction_radius * self.properties.interaction_radius;

        for id in 0..particles_map.len() as u32 {
            let mut density = 0.0;
            let mut density_near = 0.0;
            let mut neighbours = vec![];
            let p = particles_map.get(&id).unwrap();
            let pos = p.pos;
            let mut neighbours_filtered: HashMap<u32, (Vec2, f32)> = HashMap::new();

            self.get_neighbour_cell_particles(&self.grid, &p, &mut neighbours);
            for neighbour_id in neighbours.iter() {
                if *neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();

                let d = neighbour.pos - pos;

                if d.length_squared() < interaction_radius_squared {
                    let q = d.length() / self.properties.interaction_radius;

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
        dt: f32,
    ) {
        let interaction_radius_squared =
            self.properties.interaction_radius * self.properties.interaction_radius;

        for id in 0..particles_map.len() as u32 {
            let mut neighbours = vec![];
            let p = particles_map.get(&id).unwrap();
            let p_pos = p.pos;
            let p_velocity = p.velocity;

            self.get_neighbour_cell_particles(&self.grid, &p, &mut neighbours);
            for neighbour_id in neighbours.iter() {
                if *neighbour_id == id {
                    continue;
                }
                let neighbour = particles_map.get(neighbour_id).unwrap();
                let d = neighbour.pos - p_pos;

                if d.length_squared() < interaction_radius_squared {
                    let d = d.normalize();
                    let v_diff = p_velocity - neighbour.velocity;
                    let u = v_diff.dot(d);

                    if u > 0.0 {
                        let q = d.length() / self.properties.interaction_radius;
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

pub fn simulation_step(
    dt: f32,
    world: &mut ParticleWorld,
    particles_map: &mut HashMap<u32, Particle>,
    circle_shape: &CircleShape,
) {
    world.grid.neighbour_search(particles_map);
    world.viscosity(particles_map, dt);
    world.apply_gravity(particles_map, dt);
    world.predict_positions(particles_map, dt);
    world.double_density_relaxiation(particles_map, dt);
    world.check_shape_collision(particles_map, circle_shape, dt);
    world.compute_velocity_and_check_bounds(particles_map, dt);
}
