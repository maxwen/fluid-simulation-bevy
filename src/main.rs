//! Renders a 2D scene containing a single, moving sprite.

use bevy::input::mouse::MouseButtonInput;
use bevy::prelude::*;
use bevy::window::{PrimaryWindow, WindowResized, WindowResolution};
use colorgrad::{Color as GradColor, Gradient, GradientBuilder, LinearGradient};
use fluid_simulation_bevy::particle::{Particle, ParticleHashGrid, ParticleWorld};
use std::collections::HashMap;
use bevy::input::common_conditions::input_toggle_active;
use bevy_inspector_egui::quick::WorldInspectorPlugin;

const PARTICLE_AMOUNT: f32 = 2000.0;
const PARTICLE_RADIUS: f32 = 4.0;

#[derive(Copy, Clone, Component)]
struct ParticleState {
    id: u32,
}

#[derive(Component)]
struct ParticleProperties {
    color_grad: LinearGradient,
}

impl ParticleProperties {
    fn new() -> Self {
        let grad = GradientBuilder::new()
            .colors(&[
                GradColor::from((0.0, 0.0, 1.0)),
                GradColor::from((1.0, 0.0, 0.0)),
            ])
            .build::<LinearGradient>()
            .unwrap();
        ParticleProperties { color_grad: grad }
    }

    fn get_color_for_velocity(&self, velocity: f32) -> Color {
        let mapped = (1.0 / 20.0) * (f32::min(velocity, 20.0) - 0.0);
        let rgba = self.color_grad.at(mapped).to_rgba8();
        Color::srgba_u8(rgba[0], rgba[1], rgba[2], rgba[3])
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
        let world = ParticleWorld::new(width, height);
        SimulationWorld { world }
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
                resolution: WindowResolution::new(800., 600.).with_scale_factor_override(1.0),
                ..default()
            }),
            ..default()
        }))
        .add_systems(Startup, setup)
        .add_systems(FixedUpdate, update_simulation)
        .add_systems(
            Update,
            (
                update_particles,
                on_mouse_event,
                on_keyboard_event,
                on_resize_event,
            ),
        )
        .insert_resource(Time::<Fixed>::from_seconds(0.03125))
        .add_plugins(
            WorldInspectorPlugin::default().run_if(input_toggle_active(false, KeyCode::KeyI)),
        )
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
    q_window: Query<&Window, With<PrimaryWindow>>,
) {
    commands.spawn(Camera2d);

    let window = q_window.single();
    commands.spawn(SimulationGrid::new());

    let mut particles = SimulationParticles::new();
    let world = SimulationWorld::new(window.width(), window.height());
    world.world.create_particles(
        &mut particles.particles_map,
        PARTICLE_AMOUNT,
        PARTICLE_RADIUS,
    );
    let real_particle_amount = particles.particles_map.len() as f32;
    commands.spawn(particles);
    commands.spawn(world);
    commands.spawn(ParticleProperties::new());

    for id in 0..real_particle_amount as u32 {
        commands.spawn((
            Mesh2d(meshes.add(Circle::new(PARTICLE_RADIUS).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(0, 0, 255))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id },
        ));
    }
}

fn on_resize_event(mut resize_reader: EventReader<WindowResized>) {
    for e in resize_reader.read() {
        // When resolution is being changed
        println!("{:.1} x {:.1}", e.width, e.height);
    }
}

fn on_mouse_event(
    q_windows: Query<&Window, With<PrimaryWindow>>,
    mut mousebtn_evr: EventReader<MouseButtonInput>,
    grid: Query<&SimulationGrid>,
    world: Query<&SimulationWorld>,
    simulation_particles: Query<&SimulationParticles>,
) {
    use bevy::input::ButtonState;

    for ev in mousebtn_evr.read() {
        match ev.state {
            ButtonState::Pressed => {
                if let Some(position) = q_windows.single().cursor_position() {
                    let simulation_grid = grid.get_single().unwrap();
                    let simulation_world = world.get_single().unwrap();
                    let simulation_particles = simulation_particles.get_single().unwrap();

                    let mut id_list = simulation_grid
                        .particle_grid
                        .get_cell_particle_ids_from_pos(&position);
                    if id_list.len() != 0 {
                        // get any particle
                        let id = *id_list.iter().next().unwrap();
                        id_list.clear();
                        simulation_world
                            .world
                            .get_neighbour_cell_particle_ids_for_pos(
                                &simulation_particles.particles_map,
                                &simulation_grid.particle_grid,
                                id,
                                &mut id_list,
                                &position,
                            );
                        println!("{:?}", id_list);
                    }
                }
            }
            ButtonState::Released => {}
        }
    }
}

fn on_keyboard_event(
    keys: Res<ButtonInput<KeyCode>>,
    mut exit: EventWriter<AppExit>,
    mut world: Query<&mut SimulationWorld>,
    mut simulation_particles: Query<&mut SimulationParticles>,
) {
    if keys.just_pressed(KeyCode::Escape) {
        exit.send(AppExit::Success);
    }
    if keys.just_released(KeyCode::KeyQ) {
        exit.send(AppExit::Success);
    }
    if keys.just_released(KeyCode::KeyR) {
        let simulation_world = world.get_single().unwrap();
        let mut simulation_particles = simulation_particles.get_single_mut().unwrap();
        simulation_world
            .world
            .reset_particles(&mut simulation_particles.particles_map, PARTICLE_RADIUS);
    }
    if keys.just_released(KeyCode::KeyP) {
        let mut simulation_world = world.get_single_mut().unwrap();
        let mut properties = simulation_world.world.properties;
        properties.interaction_radius += 1.0;
        simulation_world
            .world
            .change_simulation_properties(properties);
    }
}

fn update_particles(
    _time: Res<Time>,
    q_window: Query<&Window, With<PrimaryWindow>>,
    camera: Query<&Transform, (With<Camera>, Without<ParticleState>)>,
    mut ui_particles: Query<(
        &mut Transform,
        &ParticleState,
        &mut MeshMaterial2d<ColorMaterial>,
    )>,
    simulation_particles: Query<&SimulationParticles>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    particle_properties: Query<&ParticleProperties>,
) {
    let window = q_window.single();
    let t = camera.get_single().unwrap();
    let particle_properties = particle_properties.get_single().unwrap();
    let simulation_particles = simulation_particles.get_single().unwrap();

    for (mut transform, p, color_handle) in &mut ui_particles {
        let simulation_particle = simulation_particles.particles_map.get(&p.id).unwrap();

        let material = materials.get_mut(&*color_handle).unwrap();
        let velocity = simulation_particle.velocity;
        material.color = particle_properties.get_color_for_velocity(velocity.length());

        let pos = particle_to_world(window, t, simulation_particle.pos);
        transform.translation.x = pos.x;
        transform.translation.y = pos.y;
    }
}

fn update_simulation(
    _fixed_time: Res<Time<Fixed>>,
    mut world: Query<&mut SimulationWorld>,
    mut grid: Query<&mut SimulationGrid>,
    mut simulation_particles: Query<&mut SimulationParticles>,
) {
    let simulation_world = world.get_single_mut().unwrap();
    let mut simulation_grid = grid.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();

    let dt = 0.3; //fixed_time.delta_secs();
    simulation_grid
        .particle_grid
        .neighbour_search(&simulation_particles.particles_map);
    simulation_world.world.viscosity(
        &mut simulation_particles.particles_map,
        &simulation_grid.particle_grid,
        dt,
    );
    simulation_world
        .world
        .apply_gravity(&mut simulation_particles.particles_map, dt);
    simulation_world
        .world
        .predict_positions(&mut simulation_particles.particles_map, dt);
    simulation_world.world.double_density_relaxiation(
        &mut simulation_particles.particles_map,
        &simulation_grid.particle_grid,
        dt,
    );
    simulation_world
        .world
        .check_boundaries_fluid(&mut simulation_particles.particles_map, dt);
    simulation_world
        .world
        .compute_velocity(&mut simulation_particles.particles_map, dt);
}
