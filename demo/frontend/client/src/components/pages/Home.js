import React, { useState } from "react";

import "../../utilities.css";
import "./Skeleton.css";

const backend_url = "http://3.101.23.89:3000/";

// Minimal webapp to send zk inputs to server generate_proof endpoint. Returns a hash. 
// Submitting the hash to the server result endpoint will return the zk proof.

const Home = () => {

  const [input, setInput] = useState("");
  const [id, setId] = useState("");
  const [proof, setProof] = useState("");

  const handleInputChange = (event) => {
    setInput(event.target.value);
  };

  const submitInput = (event) => {
    if (input !== "") {
      fetch(backend_url + "generate_proof", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: "{\"in\":" + input + "}",
      })
        .then((response) => response.json())
        .then((data) => {
          setId(data.id);
        });
    }
  }

  const submitHash = (event) => {
    if (id !== "") {
      fetch(backend_url + "result", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: "{\"id\":\"" + id + "\"}",
      })
        .then((response) => response.json())
        .then((data) => {
          setProof(JSON.stringify(data));
        });
    }
  }

  return (
    <>
    <div>
      <h1>Zk Stuff</h1>
      <div>
        <input placeholder="BinSub circuit input" type="text" value={input} onChange={handleInputChange}/>
        <button onClick={submitInput}>Submit Circuit Input</button>
      </div>
      <div><button onClick={submitHash}>Submit Hash</button></div>
      <div>{proof === "" ? "": proof}</div>
    </div>
    </>
  );
};

export default Home;
