import React from "react";
import { Router } from "@reach/router";
import NotFound from "./pages/NotFound.js";
import Skeleton from "./pages/Skeleton.js";

import "../utilities.css";

/**
 * Define the "App" component
 */
const App = () => {

  return (
    <>
      <Router>
        <Skeleton path="/" />
        <NotFound default />
      </Router>
    </>
  );
};

export default App;
