// Function to create the PMNS matrix
function createPMNSMatrix(theta12, theta23, theta13, delta) {
    const c12 = math.cos(theta12);
    const s12 = math.sin(theta12);
    const c23 = math.cos(theta23);
    const s23 = math.sin(theta23);
    const c13 = math.cos(theta13);
    const s13 = math.sin(theta13);
  
    const expNegIDelta = math.exp(math.complex(0, -delta));
    console.log(expNegIDelta);
  
    const U = math.matrix([
      [
        c12 * c13,
        s12 * c13,
        math.multiply(s13, expNegIDelta),
      ],
      [
        math.add(
          math.multiply(-s12, c23),
          math.multiply(-c12, s23, s13, expNegIDelta)
        ),
        math.add(
          math.multiply(c12, c23),
          math.multiply(-s12, s23, s13, expNegIDelta)
        ),
        math.multiply(s23, c13),
      ],
      [
        math.add(
          math.multiply(s12, s23),
          math.multiply(-c12, c23, s13, expNegIDelta)
        ),
        math.add(
          math.multiply(-c12, s23),
          math.multiply(-s12, c23, s13, expNegIDelta)
        ),
        math.multiply(c23, c13),
      ],
    ]);
  
    // Log the matrix and check if any element is NaN
    console.log('PMNS Matrix U:', U);
  
    U.forEach((val, i, j) => {
      if (math.isNaN(val)) {
        console.error(`NaN detected at U[${i}, ${j}]`);
      }
    });
  
    return U;
  }
  
  // Function to compute the transition probability
  function computeTransitionProbability(U, mSquared, L, E, alpha, beta) {
    let sum = math.complex(0, 0);
  
    // Convert energy E from MeV to GeV
    E = E / 1000;
  
    console.log('Mass squared differences:', mSquared);
    console.log('Energy (GeV):', E);
    console.log('Distance (km):', L);
  
    for (let j = 0; j < mSquared.length; j++) {
      const U_alpha_j = U.get([alpha, j]);
      const U_beta_j = U.get([beta, j]);
      const phase = math.exp(math.complex(0, (-1.267 * mSquared[j] * L) / E)); // Factor 1.267 for units: eV², km, GeV
  
      console.log(`U_alpha[${alpha}, ${j}]:`, U_alpha_j);
      console.log(`U_beta[${beta}, ${j}]:`, U_beta_j);
      console.log(`Phase for m^2[${j}] (${mSquared[j]}), L (${L}), E (${E}):`, phase);
  
      // Make sure U_alpha_j and U_beta_j are not NaN
      if (math.isNaN(U_alpha_j) || math.isNaN(U_beta_j)) {
        console.error(`NaN detected in U_alpha[${alpha}, ${j}] or U_beta[${beta}, ${j}]`);
        return NaN;
      }
  
      sum = math.add(sum, math.multiply(math.conj(U_alpha_j), U_beta_j, phase));
      console.log('Current sum:', sum);
    }
  
    console.log('Final sum:', sum);
    return math.pow(math.abs(sum), 2);
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

  function plotTransitions() {
    const theta12 = math.unit(parseFloat(document.getElementById('theta12').value), 'deg').toNumber('rad');
    const theta23 = math.unit(parseFloat(document.getElementById('theta23').value), 'deg').toNumber('rad');
    const theta13 = math.unit(parseFloat(document.getElementById('theta13').value), 'deg').toNumber('rad');
    let delta = math.unit(parseFloat(document.getElementById('delta').value), 'deg').toNumber('rad');

    const m21 = parseFloat(document.getElementById('m21').value);
    const m31 = parseFloat(document.getElementById('m31').value);
    const mSquared = [0, m21, m31];

    const E = parseFloat(document.getElementById('energy').value);

    let alpha = parseInt(document.getElementById('alpha').value);
    let isAntineutrino = false;

    // Check if the selected flavor is an antineutrino
    if (alpha >= 3) {
        alpha -= 3;
        isAntineutrino = true;
    }

    const distanceStart = parseFloat(document.getElementById('distance-start').value);
    const distanceEnd = parseFloat(document.getElementById('distance-end').value);
    const distanceStep = parseFloat(document.getElementById('distance-step').value);

    if (isAntineutrino) {
        delta = -delta;
    }

    const distances = [];
    const probabilities = [[], [], []];

    try {
        const U = createPMNSMatrix(theta12, theta23, theta13, delta);

        for (let L = distanceStart; L <= distanceEnd; L += distanceStep) {
            distances.push(L);
            for (let beta = 0; beta < 3; beta++) {
                const probability = computeTransitionProbability(U, mSquared, L, E, alpha, beta);
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
                        borderColor: 'rgba(255, 99, 132, 1)',
                        backgroundColor: 'rgba(255, 99, 132, 0.2)',
                        fill: true
                    },
                    {
                        label: labels[1],
                        data: probabilities[1],
                        borderColor: 'rgba(54, 162, 235, 1)',
                        backgroundColor: 'rgba(54, 162, 235, 0.2)',
                        fill: true
                    },
                    {
                        label: labels[2],
                        data: probabilities[2],
                        borderColor: 'rgba(75, 192, 192, 1)',
                        backgroundColor: 'rgba(75, 192, 192, 0.2)',
                        fill: true
                    }
                ]
            },
            options: {
                maintainAspectRatio: false,
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Distance (km)'
                        }
                    },
                    y: {
                        beginAtZero: true,
                        max: 1,
                        title: {
                            display: true,
                            text: 'Probability'
                        }
                    }
                },
                plugins: {
                    tooltip: {
                        mode: 'index',
                        intersect: false
                    },
                    hover: {
                        mode: 'nearest',
                        intersect: true
                    }
                }
            }
        });
    } catch (error) {
        console.error(error);
    }
    adjustCanvasSize();
  }

  document.getElementById('plot').addEventListener('click', plotTransitions);

  document.addEventListener('DOMContentLoaded', () => {
    // Adjust canvas size on initial load
    adjustCanvasSize();
    // Plot transitions on initial load
    plotTransitions();
  });

  window.addEventListener('resize', adjustCanvasSize);

