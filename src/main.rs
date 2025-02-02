//! Renders a 2D scene containing a single, moving sprite.

use bevy::ecs::system::SystemId;
use bevy::input::common_conditions::input_toggle_active;
use bevy::input::mouse::MouseButtonInput;
use bevy::prelude::*;
use bevy::window::{PrimaryWindow, WindowResized, WindowResolution};
use bevy_egui::EguiContexts;
use bevy_egui::EguiPlugin;
use colorgrad::{Color as GradColor, Gradient, GradientBuilder, LinearGradient};
use egui::Slider;
use fluid_simulation_bevy::particle::{
    simulation_step, Particle, ParticleHashGrid, ParticleWorld, ParticleWorldProperties,
};
use std::collections::HashMap;

#[derive(Copy, Clone, States, Debug, Hash, Eq, PartialEq)]
enum SimulationState {
    Running,
    Paused,
    Stepping,
}

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

#[derive(Component)]
struct SimulationProperties {
    properties: ParticleWorldProperties,
}

#[derive(Resource)]
struct CustomSystems(HashMap<String, SystemId>);

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

impl SimulationProperties {
    fn new(properties: ParticleWorldProperties) -> Self {
        SimulationProperties { properties }
    }
}
fn main() {
    let mut app = App::new();

    app.add_plugins(DefaultPlugins.set(WindowPlugin {
        primary_window: Some(Window {
            resolution: WindowResolution::new(800., 600.).with_scale_factor_override(1.0),
            ..default()
        }),
        ..default()
    }))
    .add_systems(Startup, (setup, create_particles).chain())
    .add_systems(FixedUpdate, update_simulation)
    .add_systems(
        Update,
        (
            update_particles,
            on_mouse_event,
            on_keyboard_event,
            on_resize_event,
            edit_simulation_properties.run_if(input_toggle_active(false, KeyCode::KeyP)),
        ),
    )
    .insert_resource(Time::<Fixed>::from_seconds(0.03125))
    .add_plugins(EguiPlugin)
    .insert_state(SimulationState::Running);

    let mut custom_systems = CustomSystems(HashMap::new());
    custom_systems
        .0
        .insert("reset".into(), app.register_system(recreate_particles));
    app.insert_resource(custom_systems);

    app.run();
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

    world.world.create_particles(&mut particles.particles_map);
    commands.spawn(particles);
    commands.spawn(SimulationProperties::new(world.world.properties));
    commands.spawn(world);
    commands.spawn(ParticleProperties::new());
}

fn create_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    simulation_particles: Query<&SimulationParticles>,
) {
    let simulation_particles = simulation_particles.get_single().unwrap();
    let real_particle_amount = simulation_particles.particles_map.len() as f32;

    for id in 0..real_particle_amount as u32 {
        commands.spawn((
            Mesh2d(meshes.add(Circle::new(1.0).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(0, 0, 255))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id },
        ));
    }
}

fn delete_particles(mut commands: Commands, mut particles: Query<(Entity, &Mesh2d)>) {
    for (entity, _) in &mut particles {
        commands.entity(entity).despawn();
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
                // if let Some(position) = q_windows.single().cursor_position() {
                //     let simulation_grid = grid.get_single().unwrap();
                //     let simulation_world = world.get_single().unwrap();
                //     let simulation_particles = simulation_particles.get_single().unwrap();
                //
                //     let mut id_list = simulation_grid
                //         .particle_grid
                //         .get_cell_particle_ids_from_pos(&position);
                //     if id_list.len() != 0 {
                //         // get any particle
                //         let id = *id_list.iter().next().unwrap();
                //         id_list.clear();
                //         simulation_world
                //             .world
                //             .get_neighbour_cell_particle_ids_for_pos(
                //                 &simulation_particles.particles_map,
                //                 &simulation_grid.particle_grid,
                //                 id,
                //                 &mut id_list,
                //                 &position,
                //             );
                //     }
                // }
            }
            ButtonState::Released => {}
        }
    }
}

fn on_keyboard_event(
    keys: Res<ButtonInput<KeyCode>>,
    mut exit: EventWriter<AppExit>,
    simulation_world: Query<&SimulationWorld>,
    mut simulation_particles: Query<&mut SimulationParticles>,
    state: Res<State<SimulationState>>,
    mut next_state: ResMut<NextState<SimulationState>>,
) {
    if keys.just_pressed(KeyCode::Escape) {
        exit.send(AppExit::Success);
    } else if keys.just_pressed(KeyCode::KeyQ) {
        exit.send(AppExit::Success);
    } else if keys.just_pressed(KeyCode::KeyS) {
        match state.get() {
            SimulationState::Paused | SimulationState::Running => {
                next_state.set(SimulationState::Stepping)
            }
            SimulationState::Stepping => next_state.set(SimulationState::Running),
        }
    } else if keys.just_pressed(KeyCode::Space) {
        match state.get() {
            SimulationState::Paused | SimulationState::Stepping => {
                next_state.set(SimulationState::Running)
            }
            SimulationState::Running => next_state.set(SimulationState::Paused),
        }
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
    simulation_properties: Query<&SimulationProperties>,
) {
    let window = q_window.single();
    let t = camera.get_single().unwrap();
    let particle_properties = particle_properties.get_single().unwrap();
    let simulation_particles = simulation_particles.get_single().unwrap();
    let properties = simulation_properties.get_single().unwrap();

    for (mut transform, p, color_handle) in &mut ui_particles {
        if let Some(simulation_particle) = simulation_particles.particles_map.get(&p.id) {
            let material = materials.get_mut(&*color_handle).unwrap();
            let velocity = simulation_particle.velocity;
            material.color = particle_properties.get_color_for_velocity(velocity.length());

            let pos = particle_to_world(window, t, simulation_particle.pos);
            transform.translation.x = pos.x;
            transform.translation.y = pos.y;
            transform.scale = Vec3::splat(properties.properties.particle_radius);
        }
    }
}

fn update_simulation(
    _fixed_time: Res<Time<Fixed>>,
    state: Res<State<SimulationState>>,
    mut next_state: ResMut<NextState<SimulationState>>,
    mut world: Query<&mut SimulationWorld>,
    mut grid: Query<&mut SimulationGrid>,
    mut simulation_particles: Query<&mut SimulationParticles>,
) {
    let mut simulation_world = world.get_single_mut().unwrap();
    let mut simulation_grid = grid.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();

    match state.get() {
        SimulationState::Paused => return,
        SimulationState::Running => {}
        SimulationState::Stepping => next_state.set(SimulationState::Paused),
    }

    let dt = 0.3; //fixed_time.delta_secs();
    simulation_step(
        dt,
        &mut simulation_world.world,
        &mut simulation_grid.particle_grid,
        &mut simulation_particles.particles_map,
    );
}

fn edit_simulation_properties(
    mut commands: Commands,
    mut simulation_world: Query<&mut SimulationWorld>,
    mut simulation_particles: Query<&mut SimulationParticles>,
    mut simulation_properties: Query<&mut SimulationProperties>,
    mut egui_context: EguiContexts,
    custom_systems: Res<CustomSystems>,
) {
    let mut simulation_world = simulation_world.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();
    let mut properties = simulation_properties.get_single_mut().unwrap();

    egui::Window::new("Properties").show(egui_context.ctx_mut(), |ui| {
        let mut amount = properties.properties.particle_num;
        let mut radius = properties.properties.particle_radius;
        let mut gravity = properties.properties.gravity != Vec2::ZERO;
        let mut k = properties.properties.k;
        let mut k_near = properties.properties.k_near;
        let mut rest_density = properties.properties.rest_density;
        let _interaction_radius = properties.properties.interaction_radius;
        let mut beta = properties.properties.beta;
        let mut sigma = properties.properties.sigma;

        ui.add(
            Slider::new(&mut amount, 500.0..=2000.0)
                .logarithmic(false)
                .text("Amount")
                .step_by(100.0),
        );
        properties.properties.particle_num = amount;

        ui.add(
            Slider::new(&mut radius, 2.0..=6.0)
                .logarithmic(false)
                .text("Size")
                .step_by(1.0),
        );
        properties.properties.particle_radius = radius;
        properties.properties.interaction_radius = 5.0 * radius;

        ui.checkbox(&mut gravity, "Gravity");
        if gravity {
            properties.properties.gravity = Vec2::new(0.0, 1.0);
        } else {
            properties.properties.gravity = Vec2::ZERO
        }

        ui.add(
            Slider::new(&mut rest_density, 10.0..=40.0)
                .logarithmic(false)
                .text("rest density")
                .step_by(1.0),
        );
        properties.properties.rest_density = rest_density;

        // ui.add(
        //     Slider::new(&mut interaction_radius, 10.0..=40.0)
        //         .logarithmic(false)
        //         .text("interaction radius")
        //         .step_by(1.0),
        // );
        // properties.properties.interaction_radius = interaction_radius;

        ui.add(
            Slider::new(&mut k, 0.0..=1.5)
                .logarithmic(false)
                .text("k")
                .step_by(0.01),
        );
        properties.properties.k = k;

        ui.add(
            Slider::new(&mut k_near, 0.0..=5.0)
                .logarithmic(false)
                .text("k-near")
                .step_by(0.01),
        );
        properties.properties.k_near = k_near;

        ui.add(
            Slider::new(&mut beta, 0.0..=0.5)
                .logarithmic(false)
                .text("beta")
                .step_by(0.01),
        );
        properties.properties.beta = beta;

        ui.add(
            Slider::new(&mut sigma, 0.0..=5.0)
                .logarithmic(false)
                .text("sigma")
                .step_by(0.01),
        );
        properties.properties.sigma = sigma;

        if ui.button("Apply").clicked() {
            simulation_world
                .world
                .change_simulation_properties(properties.properties);
            simulation_world
                .world
                .create_particles(&mut simulation_particles.particles_map);

            let id = custom_systems.0["reset"];
            commands.run_system(id);
        }
    });
}

fn recreate_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    simulation_particles: Query<&SimulationParticles>,
    mut particles: Query<(Entity, &Mesh2d)>,
) {
    for (entity, _) in &mut particles {
        commands.entity(entity).despawn();
    }
    let simulation_particles = simulation_particles.get_single().unwrap();
    let real_particle_amount = simulation_particles.particles_map.len() as f32;

    for id in 0..real_particle_amount as u32 {
        commands.spawn((
            Mesh2d(meshes.add(Circle::new(1.0).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(0, 0, 255))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id },
        ));
    }
}
