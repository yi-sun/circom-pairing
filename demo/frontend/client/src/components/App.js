import React from "react";
import { Router } from "@reach/router";
import NotFound from "./pages/NotFound.js";
import Home from "./pages/Home.js";

import "../utilities.css";

/**
 * Define the "App" component
 */
const App = () => {

  return (
    <>
      <Router>
        <Home path="/" />
        <NotFound default />
      </Router>
    </>
  );
};

export default App;
