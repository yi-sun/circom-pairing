import React, { useState, useEffect } from "react";

import "../../utilities.css";
import "./Skeleton.css";

// Minimal webapp to send zk inputs to server generate_proof endpoint. Returns a hash. 
// Submitting the hash to the server result endpoint will return the zk proof.

const Skeleton = () => {



  return (
    <>
      <h1>Zk Stuff</h1>
      <div>Button 1</div>
      <div>Button 2</div>
    </>
  );
};

export default Skeleton;
