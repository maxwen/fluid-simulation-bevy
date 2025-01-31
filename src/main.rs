//! Renders a 2D scene containing a single, moving sprite.

use bevy::math::bounding::{Aabb2d, Bounded2d, BoundingCircle};
use bevy::prelude::*;
use bevy::prelude::*;
use bevy::window::{PrimaryWindow, WindowResolution};
use bevyPlayground::particle::{Particle, ParticleHashGrid, ParticleWorld};
use std::collections::HashMap;

#[derive(Copy, Clone, Component)]
struct ParticleMesh {
    circle: Circle,
    radius: f32,
    id: u32,
}

#[derive(Copy, Clone, Component)]
struct ParticleState {
    id: u32,
}

impl Primitive2d for ParticleMesh {}

impl ParticleMesh {
    const fn new(id: u32, radius: f32) -> Self {
        ParticleMesh {
            circle: Circle::new(radius),
            radius,
            id,
        }
    }
}

impl Measured2d for ParticleMesh {
    fn perimeter(&self) -> f32 {
        self.circle.perimeter()
    }

    fn area(&self) -> f32 {
        self.circle.area()
    }
}

impl Bounded2d for ParticleMesh {
    fn aabb_2d(&self, isometry: impl Into<Isometry2d>) -> Aabb2d {
        self.circle.aabb_2d(isometry)
    }

    fn bounding_circle(&self, isometry: impl Into<Isometry2d>) -> BoundingCircle {
        self.circle.bounding_circle(isometry)
    }
}

impl Meshable for ParticleMesh {
    type Output = ParticleMeshBuilder;

    fn mesh(&self) -> Self::Output {
        Self::Output { particle: *self }
    }
}

struct ParticleMeshBuilder {
    particle: ParticleMesh,
}

impl MeshBuilder for ParticleMeshBuilder {
    // This is where you should build the actual mesh.
    fn build(&self) -> Mesh {
        self.particle.circle.mesh().build()
    }
}

#[derive(Component)]
struct SimulationGrid {
    particle_grid: ParticleHashGrid,
}

#[derive(Component)]
struct SimulationWorld {
    world: ParticleWorld,
}

#[derive(Component)]
struct SimulationParticles {
    particles_map: HashMap<u32, Particle>,
}

impl SimulationGrid {
    fn new() -> Self {
        SimulationGrid {
            particle_grid: ParticleHashGrid::new(),
        }
    }
}

impl SimulationWorld {
    fn new(width: f32, height: f32) -> Self {
        let world = ParticleWorld::new(width as u32, height as u32);
        SimulationWorld {
            world,
        }
    }
}

impl SimulationParticles {
    fn new() -> Self {
        SimulationParticles {
            particles_map: HashMap::new(),
        }
    }
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                resolution: WindowResolution::new(200., 200.).with_scale_factor_override(1.0),
                ..default()
            }),
            ..default()
        }))
        .add_systems(Startup, setup)
        .add_systems(FixedUpdate, update_particles)
        .run();
}

fn particle_to_world(window: &Window, camera: &Transform, position: Vec2) -> Vec2 {
    let center = camera.translation.truncate();
    let half_width = (window.width() / 2.0) * camera.scale.x;
    let half_height = (window.height() / 2.0) * camera.scale.y;
    let left = center.x - half_width;
    let top = center.y + half_height;
    Vec2::new(
        left + position.x * camera.scale.x,
        top - position.y * camera.scale.y,
    )
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    q_window: Query<&Window, With<PrimaryWindow>>
) {
    commands.spawn(Camera2d);

    let window = q_window.single();
    commands.spawn(SimulationGrid::new());

    let mut particle_amount = 500;
    let mut particles = SimulationParticles::new();
    let world  = SimulationWorld::new(window.width(), window.height());
    world.world.create_particles(&mut particles.particles_map, particle_amount);
    particle_amount = particles.particles_map.len() as u32;
    commands.spawn(particles);
    commands.spawn(world);

    for id in 0..particle_amount {
        commands.spawn((
            Mesh2d(meshes.add(ParticleMesh::new(id, 3.0).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(255, 0, 0))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id },
        ));
    }
}
fn update_particles(
    fixed_time: Res<Time<Fixed>>,
    q_window: Query<&Window, With<PrimaryWindow>>,
    mut camera: Query<&Transform, (With<Camera>, Without<ParticleState>)>,
    mut ui_particles: Query<(&mut Transform, &mut ParticleState)>,
    mut world: Query<(&mut SimulationWorld)>,
    mut grid: Query<(&mut SimulationGrid)>,
    mut simulation_particles: Query<(&mut SimulationParticles)>,

) {
    let window = q_window.single();
    let t = camera.get_single().unwrap();

    let mut simulation_world = world.get_single_mut().unwrap();
    let mut simulation_grid = grid.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();

    let dt = fixed_time.delta_secs();
    simulation_grid.particle_grid.neighbour_search(&simulation_particles.particles_map);
    simulation_world.world.apply_gravity(&mut simulation_particles.particles_map, dt);
    simulation_world.world.predict_positions(&mut simulation_particles.particles_map, dt);
    simulation_world.world.double_density_relaxiation(&mut simulation_particles.particles_map, &simulation_grid.particle_grid, dt);
    simulation_world.world.check_boundaries_fluid(&mut simulation_particles.particles_map, dt);
    simulation_world.world.compute_velocity(&mut simulation_particles.particles_map, dt);

    // let pos1 = particle_to_world(window, t, (0.0, 0.0));
    // let pos2 = particle_to_world(window, t, (100.0, 100.0));
    // let pos3 = particle_to_world(window, t, (400.0, 200.0));

    for (mut transform, p) in &mut ui_particles {
        let simulation_particle = simulation_particles.particles_map.get(&p.id).unwrap();
        let pos = particle_to_world(window, t, simulation_particle.pos);

        transform.translation.x = pos.x;
        transform.translation.y = pos.y;
    }
}
