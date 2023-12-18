import math
import time
import krpc
import numpy as np
import matplotlib.pyplot as plt


def orbital_entry(vessel, conn):
    turn_start_altitude = 250
    turn_end_altitude = 45000
    target_altitude = 190000

    # Set up streams for telemetry
    ut = conn.add_stream(getattr, conn.space_center, "ut")
    altitude = conn.add_stream(getattr, vessel.flight(), "mean_altitude")
    apoapsis = conn.add_stream(getattr, vessel.orbit, "apoapsis_altitude")
    periapsis = conn.add_stream(getattr, vessel.orbit, "periapsis_altitude")
    stage_2_resources = vessel.resources_in_decouple_stage(stage=2, cumulative=False)
    srb_fuel = conn.add_stream(stage_2_resources.amount, "SolidFuel")

    # Pre-launch setup
    vessel.control.sas = False
    vessel.control.rcs = False
    vessel.control.throttle = 1.0

    # Countdown...
    print("3...")
    time.sleep(1)
    print("2...")
    time.sleep(1)
    print("1...")
    time.sleep(1)
    print("Launch!")

    # Activate the first stage
    vessel.control.activate_next_stage()
    vessel.auto_pilot.engage()
    vessel.auto_pilot.target_pitch_and_heading(90, 90)

    # Main ascent loop
    srbs_separated = False
    turn_angle = 0
    while True:
        # Gravity turn
        if altitude() > turn_start_altitude and altitude() < turn_end_altitude:
            frac = (altitude() - turn_start_altitude) / (
                turn_end_altitude - turn_start_altitude
            )
            new_turn_angle = frac * 90
            if abs(new_turn_angle - turn_angle) > 0.5:
                turn_angle = new_turn_angle
                vessel.auto_pilot.target_pitch_and_heading(90 - turn_angle, 90)

        if vessel.thrust <= 0.000001:
            print("Stage 1 complete!")
            vessel.control.activate_next_stage()
        # Decrease throttle when approaching target apoapsis
        if apoapsis() > target_altitude * 0.9:
            vessel.auto_pilot.target_pitch_and_heading(0, 90)
            print("Approaching target apoapsis")
            break

    # Disable engines when target apoapsis is reached
    vessel.control.throttle = 0.25
    while apoapsis() < target_altitude:
        if vessel.thrust <= 0.000001:
            print("Stage 1 complete!")
            vessel.control.activate_next_stage()
        pass
    print("Target apoapsis reached")
    vessel.control.throttle = 0.0

    time.sleep(1)
    vessel.control.activate_next_stage()
    vessel.control.activate_next_stage()

    # Wait until out of atmosphere
    print("Coasting out of atmosphere")
    while altitude() < 70500:
        pass
    # Plan circularization burn (using vis-viva equation)

    print("Planning circularization burn")
    mu = vessel.orbit.body.gravitational_parameter
    r = vessel.orbit.apoapsis
    a1 = vessel.orbit.semi_major_axis
    a2 = r
    v1 = math.sqrt(mu * ((2.0 / r) - (1.0 / a1)))
    v2 = math.sqrt(mu * ((2.0 / r) - (1.0 / a2)))
    delta_v = v2 - v1
    node = vessel.control.add_node(
        ut() + vessel.orbit.time_to_apoapsis, prograde=delta_v
    )

    # Calculate burn time (using rocket equation)
    F = vessel.available_thrust
    Isp = vessel.specific_impulse * 9.82
    m0 = vessel.mass
    m1 = m0 / math.exp(delta_v / Isp)
    flow_rate = F / Isp
    burn_time = (m0 - m1) / flow_rate

    # Orientate ship
    print("Orientating ship for circularization burn")
    vessel.auto_pilot.reference_frame = node.reference_frame
    vessel.auto_pilot.target_direction = (0, 1, 0)
    vessel.auto_pilot.wait()

    # Wait until burn
    print("Waiting until circularization burn")
    burn_ut = ut() + vessel.orbit.time_to_apoapsis - (burn_time / 2.0)
    print(burn_ut)
    lead_time = 5
    conn.space_center.warp_to(burn_ut - lead_time)
    print("Warping to " + str(burn_ut - lead_time))

    # Execute burn
    print("Ready to execute burn")
    time_to_apoapsis = conn.add_stream(getattr, vessel.orbit, "time_to_apoapsis")
    while time_to_apoapsis() - (burn_time / 2.0) > 0:
        pass
    print("Executing burn")
    vessel.control.throttle = 1.0
    time.sleep(burn_time)
    vessel.control.throttle = 0.0
    vessel.control.throttle = 0.25
    while periapsis() < target_altitude:
        pass
    vessel.control.throttle = 0.0
    node.remove()

    print("Launch complete")

#doesn't work because I couldn't calculate the starting point of the maneuver
def Homan_transition(vessel, conn):
    ut = conn.add_stream(getattr, conn.space_center, "ut")
    altitude = conn.add_stream(getattr, vessel.flight(), "mean_altitude")
    apoapsis = conn.add_stream(getattr, vessel.orbit, "apoapsis_altitude")
    stage_2_resources = vessel.resources_in_decouple_stage(stage=2, cumulative=False)
    srb_fuel = conn.add_stream(stage_2_resources.amount, "SolidFuel")

    vessel.auto_pilot.disengage()

    # time.sleep(2 * math.pi * vessel.orbit.apoapsis / vessel.orbit.speed)
    print(2 * math.pi * vessel.orbit.apoapsis / vessel.orbit.speed)

    vessel = conn.space_center.active_vessel
    kerbin = conn.space_center.bodies["Kerbin"]
    mun = conn.space_center.bodies["Mun"]

    # Радиусы орбит Кербина и Муны
    r1 = kerbin.equatorial_radius + 190000  # 190 км выше поверхности Кербина
    r2 = mun.orbit.apoapsis

    # Скорости на орбитах Кербина и Муны
    v1 = vessel.orbit.speed
    v2 = mun.orbit.speed

    # Расчет delta-v для гомановского перехода
    delta_v1 = v1 * (pow((2 * r2) / (r1 + r2), 0.5) - 1)
    delta_v2 = v2 * (1 - pow((2 * r2) / (r1 + r2), 0.5))
    vessel.auto_pilot.engage()
    # Создание нода маневра
    node = vessel.control.add_node(
        ut() + vessel.orbit.time_to_apoapsis - (0.25 * vessel.orbit.period),
        prograde=delta_v1,
    )

    place = ut() + vessel.orbit.time_to_apoapsis - (0.25 * vessel.orbit.period)
    burn_ut = ut() + vessel.orbit.time_to_apoapsis + (0.25 * vessel.orbit.period)
    # Выполнение маневра
    F = vessel.available_thrust
    Isp = vessel.specific_impulse
    m0 = vessel.mass
    m1 = m0 / math.exp(delta_v1 / Isp)
    flow_rate = F / Isp
    burn_time = (m0 - m1) / flow_rate

    vessel.auto_pilot.reference_frame = node.reference_frame
    vessel.auto_pilot.target_direction = (0, 1, 0)
    vessel.auto_pilot.wait()

    print("Waiting until circularization burn")

    print(burn_ut)
    lead_time = 5
    conn.space_center.warp_to(burn_ut - (burn_time / 2.0) - lead_time)
    print("Warping to " + str(burn_ut - (burn_time / 2.0) - lead_time))
    # Выполнение маневра
    while ut() + (burn_time / 2.0) - place > 00000.1:
        pass
    print("Executing burn")
    vessel.control.throttle = 1.0
    time.sleep(burn_time)
    vessel.control.throttle = 0.0
    node.remove()

    # Удаление нода маневра
    vessel.control.remove_nodes()

    node = vessel.control.add_node(
        conn.space_center.ut + vessel.orbit.time_to_apoapsis, prograde=delta_v2
    )

    # Выполнение маневра
    vessel.auto_pilot.engage()
    vessel.auto_pilot.target_direction = node.remaining_burn_vector(
        node.reference_frame
    )
    vessel.auto_pilot.wait()
    vessel.control.throttle = 1.0
    vessel.auto_pilot.wait()
    vessel.control.throttle = 0.0

    # Удаление нода маневра
    vessel.control.remove_nodes()


def orbital_entry_mun(vessel, conn):
    turn_start_altitude = 250
    turn_end_altitude = 45000
    target_altitude = 90000

    # Set up streams for telemetry
    ut = conn.add_stream(getattr, conn.space_center, "ut")
    altitude = conn.add_stream(getattr, vessel.flight(), "mean_altitude")
    apoapsis = conn.add_stream(getattr, vessel.orbit, "apoapsis_altitude")
    periapsis = conn.add_stream(getattr, vessel.orbit, "periapsis_altitude")
    stage_2_resources = vessel.resources_in_decouple_stage(stage=2, cumulative=False)
    srb_fuel = conn.add_stream(stage_2_resources.amount, "SolidFuel")

    # Pre-launch setup
    vessel.control.sas = False
    vessel.control.rcs = False
    vessel.control.throttle = 1.0

    # Countdown...
    print("3...")
    time.sleep(1)
    print("2...")
    time.sleep(1)
    print("1...")
    time.sleep(1)
    print("Launch!")

    # Activate the first stage
    vessel.control.activate_next_stage()
    vessel.auto_pilot.engage()
    vessel.auto_pilot.target_pitch_and_heading(90, 90)

    # Main ascent loop
    srbs_separated = False
    turn_angle = 0
    while True:
        # Gravity turn
        if altitude() > turn_start_altitude and altitude() < turn_end_altitude:
            frac = (altitude() - turn_start_altitude) / (
                turn_end_altitude - turn_start_altitude
            )
            new_turn_angle = frac * 90 * (-1)
            if abs(new_turn_angle - turn_angle) > 0.5:
                turn_angle = new_turn_angle
                vessel.auto_pilot.target_pitch_and_heading(90 - turn_angle, 90)

        if vessel.thrust <= 0.000001:
            print("Stage 1 complete!")
            vessel.control.activate_next_stage()
        # Decrease throttle when approaching target apoapsis
        if apoapsis() > target_altitude * 0.9:
            vessel.auto_pilot.target_pitch_and_heading(0, 90)
            print("Approaching target apoapsis")
            break

    # Disable engines when target apoapsis is reached
    vessel.control.throttle = 0.25
    while apoapsis() < target_altitude:
        if vessel.thrust <= 0.000001:
            print("Stage 1 complete!")
            vessel.control.activate_next_stage()
        pass
    print("Target apoapsis reached")
    vessel.control.throttle = 0.0

    time.sleep(1)
    vessel.control.activate_next_stage()
    vessel.control.activate_next_stage()

    # Wait until out of atmosphere
    print("Coasting out of atmosphere")
    while altitude() < 70500:
        pass
    # Plan circularization burn (using vis-viva equation)

    print("Planning circularization burn")
    mu = vessel.orbit.body.gravitational_parameter
    r = vessel.orbit.apoapsis
    a1 = vessel.orbit.semi_major_axis
    a2 = r
    v1 = math.sqrt(mu * ((2.0 / r) - (1.0 / a1)))
    v2 = math.sqrt(mu * ((2.0 / r) - (1.0 / a2)))
    delta_v = v2 - v1
    node = vessel.control.add_node(
        ut() + vessel.orbit.time_to_apoapsis, prograde=delta_v
    )

    # Calculate burn time (using rocket equation)
    F = vessel.available_thrust
    Isp = vessel.specific_impulse * 9.82
    m0 = vessel.mass
    m1 = m0 / math.exp(delta_v / Isp)
    flow_rate = F / Isp
    burn_time = (m0 - m1) / flow_rate

    # Orientate ship
    print("Orientating ship for circularization burn")
    vessel.auto_pilot.reference_frame = node.reference_frame
    vessel.auto_pilot.target_direction = (0, 1, 0)
    vessel.auto_pilot.wait()

    # Wait until burn
    print("Waiting until circularization burn")
    burn_ut = ut() + vessel.orbit.time_to_apoapsis - (burn_time / 2.0)
    print(burn_ut)
    lead_time = 5
    conn.space_center.warp_to(burn_ut - lead_time)
    print("Warping to " + str(burn_ut - lead_time))

    # Execute burn
    print("Ready to execute burn")
    time_to_apoapsis = conn.add_stream(getattr, vessel.orbit, "time_to_apoapsis")
    while time_to_apoapsis() - (burn_time / 2.0) > 0:
        pass
    print("Executing burn")
    vessel.control.throttle = 1.0
    time.sleep(burn_time)
    vessel.control.throttle = 0.0
    vessel.control.throttle = 0.25
    while periapsis() < target_altitude:
        pass
    vessel.control.throttle = 0.0
    node.remove()

    print("Launch complete")

conn = krpc.connect(name="Launch into orbit1")
vessel = conn.space_center.active_vessel

orbital_entry(vessel, conn)

