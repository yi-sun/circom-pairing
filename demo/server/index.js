const express = require("express");
const hash = require("object-hash");
const fs = require("fs");
const { spawn } = require("child_process");
const cors = require("cors");

const app = express();

app.use(express.static("public"));

app.use(express.json());

app.get("/", (req, res) => res.send("Im a teapot"));

app.get("/slow", (req, res) => {
  // run somethign really slow here and see what happens
  var sum = 0;
  for (var i = 0; i < 10 ** 10; i++) {
    sum += i * (i - 1) * (i * 2);
  }
  res.send(`${sum}`);
});

var currentProcessesRunning = new Set();
var queue = [];

var outputData = {};

app.use(cors({ origin: "*" }));

const inQueue = (hash) => {
  // check if hash is in queue
  return queue.includes(hash);
};

const startNewProcess = (hash) => {
  if (!!outputData[hash]) return;
  const inputFileName = `inputs/${hash}.json`;
  const witnessFileName = `inputs/${hash}.wtns`;
  const publicName = `public_${hash}.json`;
  const proofName = `proof_${hash}.json`;
  // spawn a child process to run the proof generation
  const prover = spawn("sh", ["exec.sh", hash], {
    timeout: 60 * 60 * 1000,
  });
  if (!prover.pid) {
    res.status(500);
    return false;
  }
  currentProcessesRunning.add(hash);
  prover.stdout.on("data", (data) => {
    var res = data.toString();
    if (res.substring(0,8) === "{\"pi_a\":") {
        outputData[hash] = res;
        // delete the relevant files in inputs folder
        fs.unlinkSync(inputFileName);
        fs.unlinkSync(witnessFileName);
        fs.unlinkSync(publicName);
        fs.unlinkSync(proofName);
    } else {
        outputData[hash] = '{\"result\": \"BLS signature verification failed.\"}';
    }
    currentProcessesRunning.delete(hash);
    processQueue();
    prover.kill();
  });

  prover.stderr.on("data", (data) => {
    console.error(`stderr: ${data}`);
    outputData[hash] = '{\"result\": \"BLS signature verification failed.\"}';
    currentProcessesRunning.delete(hash);
    processQueue();
    prover.kill();
  });

  prover.on("close", (code) => {
    currentProcessesRunning.delete(hash);
    processQueue();
    console.log(`child process exited with code ${code}`);
  });
  return true;
};

const processQueue = () => {
  if (currentProcessesRunning.size >= 2) return null;

  // remove duplicates from queue
  queue = queue.filter(
    (item, index, self) => index === self.findIndex((t) => t === item)
  );

  // get top element from queue
  const hash = queue.shift();
  if (!hash) return;

  // check if already in process
  if (currentProcessesRunning.has(hash)) return;
  if (!!outputData[hash]) return;

  const status = startNewProcess(hash);
  return status;
};

app.post("/generate_proof", function (req, res) {
  const input = req.body;
  const inputHash = hash(input);
  const inputFileName = `inputs/${inputHash}.json`;
  console.log("input", input);
  console.log("inputFileName", inputFileName);

  if (
    currentProcessesRunning.has(inputHash) ||
    !!outputData[inputHash] ||
    inQueue(inputHash)
  ) {
    res.json({ id: inputHash });
    return;
  }

  // otherwise add to queue
  // write to file
  // and return id
  fs.writeFileSync(inputFileName, JSON.stringify(input));
  queue.push(inputHash);
  const status = processQueue();
  res.json({ id: inputHash });
});

app.post("/result", (req, res) => {
  processQueue();
  // console.log("outputData", outputData);
  const id = req.body["id"];
  const result = outputData[id];

  if (result) {
    try {
      console.log(result);
      JSON.parse(result);
      res.send(result);
    } catch (e) {
      console.log("error", e);
      res.status(404).send("{\"result\": \"ERROR\"}");
    }
  } else if (currentProcessesRunning.has(id) || inQueue(id)) {
    console.log("waiting for result");
    console.log("currentProcessesRunning", currentProcessesRunning);
    res.status(400).send("{\"result\": \"Process still running\"}");
  } else {
    console.log("ERROR! result not found", id, queue);
    // print outputdata keys
    console.log("outputData", outputData);
    res.status(404).send("{\"result\": \"ERROR\"}");
  }
});

const port = process.env.PORT || 3000;

app.listen(port, () =>
  console.log(`Server running on ${port}, http://localhost:${port}`)
);
