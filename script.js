// Function to create the PMNS matrix
function createPMNSMatrix(theta12, theta23, theta13, delta) {
    const c12 = math.cos(theta12);
    const s12 = math.sin(theta12);
    const c23 = math.cos(theta23);
    const s23 = math.sin(theta23);
    const c13 = math.cos(theta13);
    const s13 = math.sin(theta13);

    const expIDelta = math.exp(math.complex(0, delta));
    const expNegIDelta = math.conj(expIDelta);

    const U = math.matrix([
        [c12 * c13, s12 * c13, math.multiply(s13, expNegIDelta)],
        [
            math.add(math.multiply(-s12, c23), math.multiply(-c12, s23, s13, expIDelta)),
            math.add(math.multiply(c12, c23), math.multiply(-s12, s23, s13, expIDelta)),
            math.multiply(s23, c13),
        ],
        [
            math.add(math.multiply(s12, s23), math.multiply(-c12, c23, s13, expIDelta)),
            math.add(math.multiply(-c12, s23), math.multiply(-s12, c23, s13, expIDelta)),
            math.multiply(c23, c13),
        ],
    ]);

    U.forEach((val, i, j) => {
        if (math.isNaN(val)) {
            console.error(`NaN detected at U[${i}, ${j}]`);
        }
    });

    return U;
}

function computeTransMatrix(U, alpha, beta) {
    const U_star = math.conj(U);
    let re_transMatrix = math.zeros(3, 3);
    let im_transMatrix = math.zeros(3, 3);
    for (let k = 0; k < 3; k++) {
        for (let j = k + 1; j < 3; j++) {
            const product = math.multiply(U_star.get([alpha, j]), U.get([beta, j]), U.get([alpha, k]), U_star.get([beta, k]));
            const re_product = math.re(product);
            const im_product = math.im(product);
            re_transMatrix.set([j, k], re_product);
            im_transMatrix.set([j, k], im_product);
        }
    }
    return [math.multiply(-4, re_transMatrix), math.multiply(2, im_transMatrix)];
}

// Function to compute the transition probability
function computeTransitionProbability(re_transMatrix, im_transMatrix, osc_phase_exc_L, L, alpha, beta) {
    let sum = alpha === beta ? 1 : 0;
    for (let k = 0; k < 3; k++) {
        for (let j = k + 1; j < 3; j++) {
            const osc_phase = osc_phase_exc_L.get([j, k]) * L;
            const re_factor = math.pow(math.sin(osc_phase), 2);
            const im_factor = math.sin(2 * osc_phase);
            const re = re_factor * re_transMatrix.get([j, k]);
            const im = im_factor * im_transMatrix.get([j, k]);
            sum += (re + im);
        }
    }
    return sum;
}

// Update displayed values
document.querySelectorAll('input[type="range"]').forEach((slider) => {
    slider.addEventListener('input', (e) => {
        document.getElementById(`${e.target.id}-val`).textContent = e.target.value;
    });
});

let flavorChart;

function adjustCanvasSize() {
    const canvas = document.getElementById('flavorPlot');
    const container = canvas.parentElement;
    canvas.style.width = `${container.clientWidth}px`;
    canvas.style.height = `${container.clientHeight}px`;
    if (flavorChart) {
        flavorChart.resize();
    }
}

function updateHierarchy() {
    const hierarchy = document.getElementById('hierarchy').value;
    const matmLabel = document.getElementById('matm_label');

    if (hierarchy === 'normal') {
        matmLabel.textContent = '32';
    } else {
        matmLabel.textContent = '31';
    }
    checkOscillationParameters();
}

function checkOscillationParameters() {
    const hierarchy = document.getElementById('hierarchy').value;
    const matmInput = document.getElementById('matm');
    const plotButton = document.getElementById('plot');
    const oscillationFieldset = document.querySelector('fieldset:nth-of-type(1)'); // Assuming Oscillation Parameters is the first fieldset
    let warning = document.getElementById('matm-warning');

    if ((hierarchy === 'normal' && parseFloat(matmInput.value) < 0) || (hierarchy === 'inverted' && parseFloat(matmInput.value) > 0)) {
        if (!warning) {
            warning = document.createElement('div');
            warning.id = 'matm-warning';
            warning.style.color = 'red';
            warning.style.fontWeight = 'bold';
            warning.textContent = 'Warning: Δm_atm² has the wrong sign for the selected hierarchy.';
            oscillationFieldset.insertBefore(warning, oscillationFieldset.firstChild);
        }
        matmInput.style.borderColor = 'red';
        plotButton.disabled = true;
    } else {
        if (warning) {
            warning.remove();
        }
        matmInput.style.borderColor = '';
        plotButton.disabled = false;
    }
}

function checkPlotPoints() {
    const distanceStart = parseFloat(document.getElementById('distance-start').value);
    const distanceEnd = parseFloat(document.getElementById('distance-end').value);
    const distanceStep = parseFloat(document.getElementById('distance-step').value);
    const plotOptionsFieldset = document.querySelector('fieldset:nth-of-type(3)'); // Assuming Plot Options is the third fieldset

    const numPoints = (distanceEnd - distanceStart) / distanceStep;
    console.log("numPoints: ", numPoints);
    let warning = document.getElementById('plot-warning');

    if (numPoints > 5000) {
        if (!warning) {
            warning = document.createElement('div');
            warning.id = 'plot-warning';
            warning.style.color = 'red';
            warning.style.fontWeight = 'bold';
            warning.textContent = 'Warning: The specified settings will compute more than 5,000 points per flavour.';
            plotOptionsFieldset.insertBefore(warning, plotOptionsFieldset.firstChild);
        }
    } else {
        if (warning) {
            warning.remove();
        }
    }
}

function plotTransitions() {
    const theta12 = math.unit(parseFloat(document.getElementById('theta12').value), 'deg').toNumber('rad');
    const theta23 = math.unit(parseFloat(document.getElementById('theta23').value), 'deg').toNumber('rad');
    const theta13 = math.unit(parseFloat(document.getElementById('theta13').value), 'deg').toNumber('rad');
    let delta = math.unit(parseFloat(document.getElementById('delta').value), 'deg').toNumber('rad');

    const matm = parseFloat(document.getElementById('matm').value);
    const msol = parseFloat(document.getElementById('msol').value);
    const hierarchy = document.getElementById('hierarchy').value;
    const isStacked = document.getElementById('stacked').value === 'true';
    
    if (hierarchy === 'normal' && matm < 0) {
        alert("For normal hierarchy, Δm_atm² should be > 0.");
        return;
    } else if (hierarchy === 'inverted' && matm > 0) {
        alert("For inverted hierarchy, Δm_atm² should be < 0.");
        return;
    }

    let mSquared;
    if (hierarchy === 'normal') {
        // Normal hierarchy
        mSquared = math.matrix([
            [0, msol, matm + msol],
            [-msol, 0, matm],
            [-(matm + msol), -matm, 0]
        ]);
    } else {
        // Inverted hierarchy
        mSquared = math.matrix([
            [0, msol, matm],
            [-msol, 0, matm - msol],
            [-matm, -(matm - msol), 0]
        ]);
    }

    let E = parseFloat(document.getElementById('energy').value);
    // Convert energy E from MeV to GeV
    E /= 1000; 

    let alpha = parseInt(document.getElementById('alpha').value);
    let isAntineutrino = false;

    // Check if the selected flavor is an antineutrino
    if (alpha >= 3) {
        alpha -= 3;
        isAntineutrino = true;
        delta *= -1;
    }

    const distanceStart = parseFloat(document.getElementById('distance-start').value);
    const distanceEnd = parseFloat(document.getElementById('distance-end').value);
    const distanceStep = parseFloat(document.getElementById('distance-step').value);

    const probabilities = [[], [], []];
    const osc_phase_exc_L = math.multiply(1.27 / E, mSquared);
    const distances = Array.from(
      { length: Math.ceil((distanceEnd - distanceStart) / distanceStep) + 1 },
      (_, i) => distanceStart + i * distanceStep
    );
    const U = createPMNSMatrix(theta12, theta23, theta13, delta);

    for (let beta = 0; beta < 3; beta++) {
      const [re_transMatrix, im_transMatrix] = computeTransMatrix(U, alpha, beta);
      for (let L of distances) {
        const probability = computeTransitionProbability(re_transMatrix, im_transMatrix, osc_phase_exc_L, L, alpha, beta);
        probabilities[beta].push(probability);
      }
    }

    const ctx = document.getElementById('flavorPlot').getContext('2d');
    
    // Destroy existing chart instance if it exists
    if (flavorChart) {
        flavorChart.destroy();
    }

    const labels = isAntineutrino ? ['ν̅e', 'ν̅μ', 'ν̅τ'] : ['νe', 'νμ', 'ντ'];

    flavorChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: distances,
            datasets: [
                {
                    label: labels[0],
                    data: probabilities[0],
                    borderColor: isStacked ? 'rgba(255, 99, 132, 0)' : 'rgba(255, 99, 132, 1)',
                    backgroundColor: isStacked ? 'rgba(255, 99, 132, 0.7)' : 'rgba(255, 99, 132, 0.2)',
                    fill: true,
                    pointRadius: 0
                },
                {
                    label: labels[1],
                    data: probabilities[1],
                    borderColor: isStacked ? 'rgba(54, 162, 235, 0)' : 'rgba(54, 162, 235, 1)',
                    backgroundColor: isStacked ? 'rgba(54, 162, 235, 0.7)' : 'rgba(54, 162, 235, 0.2)',
                    fill: true,
                    pointRadius: 0
                },
                {
                    label: labels[2],
                    data: probabilities[2],
                    borderColor: isStacked ? 'rgba(75, 192, 192, 0)' : 'rgba(75, 192, 192, 1)',
                    backgroundColor: isStacked ? 'rgba(75, 192, 192, 0.7)' : 'rgba(75, 192, 192, 0.2)',
                    fill: true,
                    pointRadius: 0
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: (ctx) => 'Neutrino Transition Probabilities',
                    font: {
                        size: 24 // Double the default font size
                    }
                },
                tooltip: {
                    mode: 'index'
                },
                legend: {
                    labels: {
                        font: {
                            size: 18 // Increase legend font size
                        }
                    }
                }
            },
            interaction: {
                mode: 'nearest',
                axis: 'x',
                intersect: false
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Distance (km)',
                        font: {
                            size: 24 // Double the default font size
                        }
                    },
                    ticks: {
                        font: {
                            size: 18 // Increase axis label font size
                        }
                    }
                },
                y: {
                    stacked: isStacked,
                    title: {
                        display: true,
                        text: 'Probability',
                        font: {
                            size: 24 // Double the default font size
                        }
                    },
                    ticks: {
                        font: {
                            size: 18 // Increase axis label font size
                        }
                    },
                    beginAtZero: true,
                    max: 1
                }
            }
        }
    });
    adjustCanvasSize();
}

document.getElementById('plot').addEventListener('click', plotTransitions);
document.getElementById('hierarchy').addEventListener('change', updateHierarchy);
document.getElementById('hierarchy').addEventListener('change', checkOscillationParameters);
document.getElementById('matm').addEventListener('input', checkOscillationParameters);
document.getElementById('stacked').addEventListener('change', adjustCanvasSize);
document.getElementById('distance-start').addEventListener('input', checkPlotPoints);
document.getElementById('distance-end').addEventListener('input', checkPlotPoints);
document.getElementById('distance-step').addEventListener('input', checkPlotPoints);

document.addEventListener('DOMContentLoaded', () => {
    // Adjust canvas size on initial load
    adjustCanvasSize();
    // Plot transitions on initial load
    plotTransitions();
});

window.addEventListener('resize', adjustCanvasSize);

