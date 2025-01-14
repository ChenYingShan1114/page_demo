import {OrbitControls} from 'https://unpkg.com/three@0.127.0/examples/jsm/controls/OrbitControls.js'
import * as THREE from 'https://unpkg.com/three@0.127.0/build/three.module.js';

// constant
const eV = 1.60217662e-19; // J
const qe = 1.602176621e-19; // C
const eps0 = 8.854187817e-12; // F/m
const me = 9.1093837e-31; // kg
const pi = 3.1415926535899;
const c = 299792458;

// normalization constant
const nor_position = 1e-15;
const nor_time = nor_position / c;
const nor_velocity = c;
const nor_charge = qe;
const nor_ke = 4 * pi * eps0;
const nor_mass = me;
const nor_energy = nor_mass * nor_velocity * nor_velocity;
const nor_force = nor_charge * nor_charge / nor_ke / nor_position / nor_position;

const nor_const = nor_charge * nor_charge / nor_velocity / nor_velocity / nor_mass / nor_position / nor_ke;

// alpha particle
let q1 = 2 * qe / nor_charge;
let m1 = 6.644657230e-27 / nor_mass;
let E_in = 5 * 1e6 * eV / nor_energy;
let v_in = Math.sqrt(2 * E_in / m1);
// golden nucleus
let q2 = 79 * qe / nor_charge;
// simulation parameter
let ini_x = -1e-12 / nor_position;
let total_time = 2 * Math.abs(ini_x) / v_in; // unit: nor_time  
let steps = 300;
const ke = q1 * q2;

function ax(x, y, vx, vy, time){
    let dist = Math.sqrt(x*x+y*y);
    return ke / m1 / Math.pow(dist,3) * (x - 0) * nor_const; 
}
function ay(x, y, vx, vy, time){
    let dist = Math.sqrt(x*x+y*y);
    return ke / m1 / Math.pow(dist,3) * (y - 0) * nor_const; 
}

let x_rk4_new = 0;
let y_rk4_new = 0;
let vx_rk4_new = 0;
let vy_rk4_new = 0;
function rk4(x_old, y_old, vx_old, vy_old, time, tau){

    let  x_RK1,  x_RK2,  x_RK3,  x_RK4,  y_RK1,  y_RK2,  y_RK3,  y_RK4;
    let vx_RK1, vx_RK2, vx_RK3, vx_RK4, ax_RK1, ax_RK2, ax_RK3, ax_RK4;
    let vy_RK1, vy_RK2, vy_RK3, vy_RK4, ay_RK1, ay_RK2, ay_RK3, ay_RK4;

     x_RK1 =  x_old;
     y_RK1 =  y_old;
    vx_RK1 = vx_old;
    vy_RK1 = vy_old;
    ax_RK1 = ax(x_RK1, y_RK1, vx_RK1, vy_RK1, time);
    ay_RK1 = ay(x_RK1, y_RK1, vx_RK1, vy_RK1, time);

     x_RK2 =  x_old + 0.5 * tau * vx_RK1;
     y_RK2 =  y_old + 0.5 * tau * vy_RK1;
    vx_RK2 = vx_old + 0.5 * tau * ax_RK1;
    vy_RK2 = vy_old + 0.5 * tau * ay_RK1;
    ax_RK2 = ax(x_RK2, y_RK2, vx_RK2, vy_RK2, time + 0.5 * tau);
    ay_RK2 = ay(x_RK2, y_RK2, vx_RK2, vy_RK2, time + 0.5 * tau);

     x_RK3 =  x_old + 0.5 * tau * vx_RK2;
     y_RK3 =  y_old + 0.5 * tau * vy_RK2;
    vx_RK3 = vx_old + 0.5 * tau * ax_RK2;
    vy_RK3 = vy_old + 0.5 * tau * ay_RK2;
    ax_RK3 = ax(x_RK3, y_RK3, vx_RK3, vy_RK3, time + 0.5 * tau);
    ay_RK3 = ay(x_RK3, y_RK3, vx_RK3, vy_RK3, time + 0.5 * tau);

     x_RK4 =  x_old + tau * vx_RK3;
     y_RK4 =  y_old + tau * vy_RK3;
    vx_RK4 = vx_old + tau * ax_RK3;
    vy_RK4 = vy_old + tau * ay_RK3;
    ax_RK4 = ax(x_RK4, y_RK4, vx_RK4, vy_RK4, time + tau);
    ay_RK4 = ay(x_RK4, y_RK4, vx_RK4, vy_RK4, time + tau);

    vx_rk4_new = vx_old + tau / 6 * (ax_RK1 + 2 * ax_RK2 + 2 * ax_RK3 + ax_RK4);
    vy_rk4_new = vy_old + tau / 6 * (ay_RK1 + 2 * ay_RK2 + 2 * ay_RK3 + ay_RK4);
     x_rk4_new =  x_old + tau / 6 * (vx_RK1 + 2 * vx_RK2 + 2 * vx_RK3 + vx_RK4);
     y_rk4_new =  y_old + tau / 6 * (vy_RK1 + 2 * vy_RK2 + 2 * vy_RK3 + vy_RK4);

}

let x_rka_new = 0;
let y_rka_new = 0;
let vx_rka_new = 0;
let vy_rka_new = 0;
let tau_rka = 0;
let dt_rka = 0;
function rka(x_old, y_old, vx_old, vy_old, time, tau) {
    let vx_small_1 = 0, vy_small_1 = 0, x_small_1 = 0, y_small_1 = 0;
    let vx_small_2 = 0, vy_small_2 = 0, x_small_2 = 0, y_small_2 = 0;
    let vx_big     = 0, vy_big     = 0, x_big     = 0, y_big     = 0;
    let r_small = 0, r_big = 0;
    for (let j = 0; j <= 100; j++){ 
        // two small steps
        rk4(x_old, y_old, vx_old, vy_old, time, 0.5 * tau);
        vx_small_1 = vx_rk4_new;
        vy_small_1 = vy_rk4_new;
        x_small_1 = x_rk4_new;
        y_small_1 = y_rk4_new;
        rk4(x_small_1, y_small_1, vx_small_1, vy_small_1, time + 0.5 * tau, 0.5 * tau);
        vx_small_2 = vx_rk4_new;
        vy_small_2 = vy_rk4_new;
         x_small_2 = x_rk4_new;
         y_small_2 = y_rk4_new;
         r_small = Math.sqrt(x_small_2 * x_small_2 + y_small_2 * y_small_2);
        
        // one big steps
        rk4(x_old, y_old, vx_old, vy_old, time, tau);
        vx_big = vx_rk4_new;
        vy_big = vy_rk4_new;
         x_big = x_rk4_new;
         y_big = y_rk4_new;
         r_big = Math.sqrt(x_big * x_big + y_big * y_big);
        
        let s1 = 0.99, s2 = 1.5;
        let ideal_err = 1e-10;
        // 若local_err > ideal_err大，dt就變小。
        // 若local_err < ideal_err小，dt就變大
        let local_err = Math.abs((r_small - r_big) / r_big);

        let dt_est = tau * Math.pow(ideal_err / local_err, 0.2);
        let dt_new = s1 * dt_est;
        if (s1 * dt_est > s2 * tau){
            dt_new = s2 * tau;
        }
        else if (s1 * dt_est < tau / s2){
            dt_new = tau / s2;
        }
        else{
            dt_new = s1 * dt_est;
        }
        tau = dt_new;
        if (local_err < ideal_err){
            //vx.push(vx_big);
            //vy.push(vy_big);
            //x.push(x_big);
            //y.push(x_big);
            x_rka_new = x_big;
            y_rka_new = y_big;
            vx_rka_new = vx_big;
            vy_rka_new = vy_big;
            break;
        }            
    }
    //t_now.push(t_now[step-1] + dt);
    tau_rka = time + tau;
    dt_rka = tau;

}

let clock =  new THREE.Clock();

const canvas = document.querySelector('canvas.webgl')

// Scene
const scene = new THREE.Scene();
// Camera
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 100);

const renderer = new THREE.WebGLRenderer({
    canvas: canvas,
    alpha: true,
})
renderer.shadowMap.enabled = true;
renderer.setSize( window.innerWidth, window.innerHeight );
renderer.setAnimationLoop( animate );
//document.body.appendChild( canvas );

// Controls
const orbit = new OrbitControls(camera, canvas);
orbit.target.set(0, 0, 0);

const axesHelper = new THREE.AxesHelper( 1000 );
scene.add( axesHelper );

const sphere_au = new THREE.Mesh( new THREE.SphereGeometry( 25, 32, 16 ), new THREE.MeshStandardMaterial( { color: 0xaaaa00 } ) );
scene.add( sphere_au );
sphere_au.receiveShadow = true;

const ambientLight = new THREE.AmbientLight(0xffffff);
scene.add( ambientLight );

/*
const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
scene.add( directionalLight );
directionalLight.position.set(-50, 50, 0);
directionalLight.castShadow = true;

const dLightHelper = new THREE.DirectionalLightHelper(directionalLight, 5);
scene.add(dLightHelper);
const dLightShadowHelper = new THREE.CameraHelper(directionalLight.shadow.camera);
scene.add(dLightShadowHelper);
*/

const spotLight = new THREE.SpotLight(0xffffff);
scene.add(spotLight);
spotLight.position.set(-50, 50, 0);
spotLight.castShadow = true;
spotLight.intensity = 10000;

const sLightHelper = new THREE.SpotLightHelper(spotLight);
scene.add(sLightHelper);

let num = 200;
let vx = [], vy = [], x = [], y = [], t_now = [], dt = [], spheres = [];
let index = -1;
let dis = 500;
for (let i = 0; i < num; i++) {
    const ini_y = Math.pow(-1, i) * Math.random() * 100 * 1e-15 / nor_position;

    vx.push(v_in);
    vy.push(0);
    x.push(ini_x);
    y.push(ini_y);
    t_now.push(0);
    dt.push(3000);

    const sphere = new THREE.Mesh( new THREE.SphereGeometry( 1, 32, 16 ), new THREE.MeshStandardMaterial( { color: 0xaaaaaa } ) );
    scene.add( sphere );
    sphere.position.set(ini_x, ini_y, 0);
    spheres.push(sphere);
    sphere.castShadow = true;

    if (ini_y < dis) {
        dis = ini_y;
        index = i;
    }
}

camera.position.set(ini_x, 6, 26);
orbit.target.set(x[index], 0, 0);
orbit.update();

let step = 1;


function animate() {
    const delta = clock.getElapsedTime();

    //rka
    if (step < 25 * delta && step < steps) {

        for (let i = 0; i < num; i++) {
            rka(x[i], y[i], vx[i], vy[i], t_now[i], dt[i]);
            vx[i] = vx_rka_new; vy[i] = vy_rka_new; x[i] = x_rka_new; y[i] = y_rka_new; t_now[i] = tau_rka; dt[i] = dt_rka;
            spheres[i].position.x = x[i];
            spheres[i].position.y = y[i];
            spheres[i].material.color.setRGB( 0.3, 1 - vx[i] / v_in, 1 );

            const sphere_tmp = new THREE.Mesh( new THREE.SphereGeometry( 0.5, 8, 4 ), new THREE.MeshStandardMaterial( { color: 0x000000 } ) );
            sphere_tmp.material.color.setRGB( 0.3, 1 - vx[i] / v_in, 1 );
            //scene.add( sphere_tmp );
            sphere_tmp.position.set(x[i], y[i], 0);
        }

        step++;
        if ( x[index] < -30) {
            camera.position.x = x[index];
            orbit.target.x = x[index];
        } else {
            camera.position.z += 0.5;
        }
	    orbit.update();
    }

    /*rk4
    if (step < 50 * delta && step < steps) {
        rk4(x[step-1], y[step-1], vx[step-1], vy[step-1], t_now[step-1], dt);
        vx.push(vx_new);
        vy.push(vy_new);
        x.push(x_new);
        y.push(y_new);
        t_now.push(t_now[step-1] + dt);
        console.log(t_now[step], x[step], y[step]);
        sphere_e.position.x = x[step];
        sphere_e.position.y = y[step];
        step++;
    }*/
    
    renderer.render( scene, camera );

}
