//! Renders a 2D scene containing a single, moving sprite.

use bevy::ecs::system::SystemId;
use bevy::input::common_conditions::input_toggle_active;
use bevy::input::mouse::MouseButtonInput;
use bevy::prelude::*;
use bevy::window::{PrimaryWindow, WindowResolution};
use bevy_egui::EguiContexts;
use bevy_egui::EguiPlugin;
use bevy_tokio_tasks::{TokioTasksPlugin, TokioTasksRuntime};
use colorgrad::{Color as GradColor, Gradient, GradientBuilder, LinearGradient};
use egui::Slider;
use fluid_simulation_bevy::particle::{
    simulation_step, CircleShape, Particle, ParticleWorld, ParticleWorldProperties,
};
use std::collections::HashMap;
use std::time::Duration;

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
    properties: ParticleWorldProperties,
}

impl ParticleProperties {
    fn new() -> Self {
        let grad = GradientBuilder::new()
            .colors(&[
                GradColor::from((0.0, 0.0, 1.0)),
                GradColor::from((1.0, 1.0, 0.0)),
                GradColor::from((1.0, 0.0, 0.0)),
            ])
            .build::<LinearGradient>()
            .unwrap();
        ParticleProperties {
            color_grad: grad,
            properties: ParticleWorldProperties::new(),
        }
    }

    fn get_color_for_velocity(&self, velocity_squared: f32) -> Color {
        let mapped = (1.0 / 200.0) * (f32::min(velocity_squared, 200.0) - 0.0);
        let rgba = self.color_grad.at(mapped).to_rgba8();
        Color::srgba_u8(rgba[0], rgba[1], rgba[2], rgba[3])
    }
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
struct SimulationActivators {
    circle: CircleShape,
}

#[derive(Resource)]
struct CustomSystems(HashMap<String, SystemId>);

impl SimulationWorld {
    fn new(width: f32, height: f32, properties: ParticleWorldProperties) -> Self {
        let world = ParticleWorld::new(width, height, properties);
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

impl SimulationActivators {
    fn new() -> Self {
        SimulationActivators {
            circle: CircleShape::new(Vec2::ZERO, 40.0, false),
        }
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
    .add_systems(
        Startup,
        (
            (setup_simulation, create_particles).chain(),
            tokio_background_task,
        ),
    )
    .add_systems(FixedUpdate, update_simulation)
    .add_systems(
        Update,
        (
            update_particles,
            on_mouse_event,
            on_keyboard_event,
            edit_simulation_properties.run_if(input_toggle_active(false, KeyCode::KeyP)),
        ),
    )
    .insert_resource(Time::<Fixed>::from_seconds(0.03125))
    .add_plugins(EguiPlugin)
    .insert_state(SimulationState::Running)
    .add_plugins(TokioTasksPlugin::default());

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

fn setup_simulation(mut commands: Commands, window: Query<&Window, With<PrimaryWindow>>) {
    let particle_properties = ParticleProperties::new();
    let world = SimulationWorld::new(
        window.single().width(),
        window.single().height(),
        particle_properties.properties,
    );
    let mut particles = SimulationParticles::new();
    world.world.create_particles(&mut particles.particles_map);

    commands.spawn(particles);
    commands.spawn(world);
    commands.spawn(particle_properties);
    commands.spawn(SimulationActivators::new());
}

fn create_particles(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    simulation_particles: Query<&SimulationParticles>,
) {
    commands.spawn(Camera2d);

    let simulation_particles = simulation_particles.get_single().unwrap();
    let real_particle_amount = simulation_particles.particles_map.len();

    for id in 0..real_particle_amount {
        commands.spawn((
            Mesh2d(meshes.add(Circle::new(1.0).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(0, 0, 255))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id: id as u32 },
        ));
    }
}

#[allow(dead_code)]
fn delete_particles(mut commands: Commands, mut particles: Query<(Entity, &Mesh2d)>) {
    for (entity, _) in &mut particles {
        commands.entity(entity).despawn();
    }
}

// fn on_resize_event(mut resize_reader: EventReader<WindowResized>) {
//     for e in resize_reader.read() {
//         // When resolution is being changed
//         println!("{:.1} x {:.1}", e.width, e.height);
//     }
// }

fn on_mouse_event(
    window: Query<&Window, With<PrimaryWindow>>,
    mut mouse_events: EventReader<MouseButtonInput>,
    mut simulation_activators: Query<&mut SimulationActivators>,
    keys: Res<ButtonInput<KeyCode>>,
) {
    use bevy::input::ButtonState;
    let mut simulation_activators = simulation_activators.get_single_mut().unwrap();

    for ev in mouse_events.read() {
        match ev.state {
            ButtonState::Pressed => {
                let shift = keys.any_pressed([KeyCode::ShiftLeft, KeyCode::ShiftRight]);
                simulation_activators.circle.enable(shift);
            }
            ButtonState::Released => {
                simulation_activators.circle.disable();
            }
        }
    }

    if let Some(position) = window.single().cursor_position() {
        if simulation_activators.circle.is_active() {
            simulation_activators.circle.set_position(position);
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
    } else if keys.just_pressed(KeyCode::KeyR) {
        let simulation_world = simulation_world.get_single().unwrap();
        let mut simulation_particles = simulation_particles.get_single_mut().unwrap();
        simulation_world
            .world
            .create_particles(&mut simulation_particles.particles_map);
    }
}

fn update_particles(
    window: Query<&Window, With<PrimaryWindow>>,
    camera: Query<&Transform, (With<Camera>, Without<ParticleState>)>,
    mut ui_particles: Query<(
        &mut Transform,
        &ParticleState,
        &mut MeshMaterial2d<ColorMaterial>,
    )>,
    mut materials: ResMut<Assets<ColorMaterial>>,
    particle_properties: Query<&ParticleProperties>,
    simulation_particles: Query<&SimulationParticles>,
) {
    let window = window.single();
    let camera = camera.get_single().unwrap();
    let particle_properties = particle_properties.get_single().unwrap();
    let simulation_particles = simulation_particles.get_single().unwrap();

    let particle_radius_scale = Vec3::splat(particle_properties.properties.particle_radius);
    for (mut transform, p, color_handle) in &mut ui_particles {
        if let Some(simulation_particle) = simulation_particles.particles_map.get(&p.id) {
            let material = materials.get_mut(&*color_handle).unwrap();
            let velocity = simulation_particle.velocity;
            material.color = particle_properties.get_color_for_velocity(velocity.length_squared());

            let pos = particle_to_world(window, camera, simulation_particle.pos);
            transform.translation.x = pos.x;
            transform.translation.y = pos.y;
            transform.scale = particle_radius_scale;
        }
    }
}

fn update_simulation(
    _fixed_time: Res<Time<Fixed>>,
    state: Res<State<SimulationState>>,
    mut next_state: ResMut<NextState<SimulationState>>,
    mut world: Query<&mut SimulationWorld>,
    mut simulation_particles: Query<&mut SimulationParticles>,
    simulation_activators: Query<&SimulationActivators>,
) {
    let mut simulation_world = world.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();
    let simulation_activators = simulation_activators.get_single().unwrap();

    match state.get() {
        SimulationState::Paused => return,
        SimulationState::Running => {}
        SimulationState::Stepping => next_state.set(SimulationState::Paused),
    }

    let dt = 0.5; //fixed_time.delta_secs();
    simulation_step(
        dt,
        &mut simulation_world.world,
        &mut simulation_particles.particles_map,
        &simulation_activators.circle,
    );
}

fn tokio_background_task(runtime: ResMut<TokioTasksRuntime>) {
    runtime.spawn_background_task(|mut ctx| async move {
        loop {
            ctx.run_on_main_thread(move |ctx| {
                println!("Check task");
                if let Some(state) = ctx.world.get_resource::<State<SimulationState>>() {
                    match state.get() {
                        SimulationState::Paused => {
                            println!("Paused")
                        }
                        SimulationState::Running => {
                            println!("Running")
                        }
                        SimulationState::Stepping => {
                            println!("Stepping")
                        }
                    }
                }
            })
            .await;
            tokio::time::sleep(Duration::from_secs(1)).await;
        }
    });
}

fn edit_simulation_properties(
    mut commands: Commands,
    mut simulation_world: Query<&mut SimulationWorld>,
    mut simulation_particles: Query<&mut SimulationParticles>,
    mut particle_properties: Query<&mut ParticleProperties>,
    mut egui_context: EguiContexts,
    custom_systems: Res<CustomSystems>,
) {
    let mut simulation_world = simulation_world.get_single_mut().unwrap();
    let mut simulation_particles = simulation_particles.get_single_mut().unwrap();
    let mut properties = particle_properties.get_single_mut().unwrap();

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
            Slider::new(&mut amount, 500.0..=3000.0)
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
            Slider::new(&mut rest_density, 0.0..=40.0)
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

        // let serialized = serde_json::to_string(&properties.properties).unwrap();
        // println!("{}", serialized);
        // let deserialized: ParticleWorldProperties = serde_json::from_str(&serialized).unwrap();

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
    let real_particle_amount = simulation_particles.particles_map.len();

    for id in 0..real_particle_amount {
        commands.spawn((
            Mesh2d(meshes.add(Circle::new(1.0).mesh())),
            MeshMaterial2d(materials.add(Color::srgb_u8(0, 0, 255))),
            Transform::from_xyz(0.0, 0.0, 0.0),
            ParticleState { id: id as u32 },
        ));
    }
}
