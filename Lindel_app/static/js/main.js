function showExample() {
    var exampleSeq = 'TAACGTTATCAACGCCTATATTAAAGCGACCGTCGGTTGAACTGCGTGGATCAATGCGTC';
    document.getElementById('seq').value = exampleSeq;
}

function submit_sequence() {
    // Clear previous results and errors
    document.getElementById('error-message').style.display = 'none';
    Plotly.purge('graph');
    document.getElementById('blank-space').style.display = 'block';

    $.ajax({
        method: 'POST',
        url: "http://127.0.0.1:5000/predict",
        data: {'string': document.getElementById('seq').value},
        success: function(response) {
            obj = JSON.parse(response);
            $(document).ready(function () {
                // Hide the blank space div
                document.getElementById('blank-space').style.display = 'none';
                plottable(obj);
            });
        },
        error: function(response) {
            console.error(response);
            // Display error message
            var errorMessageDiv = document.getElementById('error-message');
            errorMessageDiv.style.display = 'block';
            var errorResponse = JSON.parse(response.responseText);
            errorMessageDiv.innerText = 'Error: ' + errorResponse.error;
        }
    });
}

var plottable = function(plotobj) {
    var header = ['Sequence', 'Frequency', 'Indels'];
    var sequence = plotobj[0][header[0]];
    var dict = {};
    dict[header[0]] = [];
    dict[header[1]] = [];
    dict[header[2]] = [];
    dict['table'] = [];
    for (var i = 1; i < 21; i++) {
        header.forEach(function(key) {
            dict[key].push(plotobj[i][key]);
            if (key == 'Indels') {
                dict['table'].push(plotobj[i][key].split(' ')[0]);
            }
        });
    }
    var values = [dict[header[0]], dict[header[1]], dict['table']];
    var table = {
        type: 'table',
        columnwidth: [600, 150, 100],
        header: {
            values: [[sequence], ['Frequency (%)'], ['Indels']],
            align: "center",
            height: 20,
            line: {width: 0, color: 'black'},
            fill: {color: "grey"},
            font: {family: "Courier New, monospace", size: '12pt', color: "white"}
        },
        cells: {
            values: values,
            align: "center",
            height: 20,
            line: {color: "black", width: 0},
            font: {family: "Courier New, monospace", size: '12pt', color: ["black"]}
        },
        domain: {x: [0, 0.65], y: [0, 1]},
    };

    var data = {
        y: dict[header[2]],
        x: dict[header[1]],
        type: 'bar',
        orientation: 'h',
        xaxis: 'x1',
        yaxis: 'y1',
        width: 0.8,
        name: '',
        marker: {
            color: '#c4c4c4'
        }
    };
    var total = [table, data];
    var layout = {
        width: 1200,
        height: 600,
        showlegend: false,
        bargap: 0.1,
        hovermode: 'closest',
        xaxis: {domain: [0.7, 1], anchor: 'y1', showticklabels: true, fixedrange: true},
        yaxis: {domain: [0, 0.95], autorange: 'reversed', anchor: 'x1', showticklabels: false, fixedrange: true},
    };
    Plotly.plot('graph', total, layout);
}