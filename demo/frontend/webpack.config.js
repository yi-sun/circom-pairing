/*
|--------------------------------------------------------------------------
| webpack.config.js -- Configuration for Webpack
|--------------------------------------------------------------------------
|
| Webpack turns all the clientside HTML, CSS, Javascript into one bundle.js file.
| This is done for performance reasons, as well as for compatability reasons.
|
| You do not have to worry about this file, except for proxy section below.
| All proxies does is route traffic from the hotloader to the backend.
| You must define explicity all routes here, as we do for the /api/* routes.
|
| The rest of this file tell webpack which types of files to bundle (in the rules).
| In addition, it also uses babel to transpile your javascript into code all browsers can use.
| see https://babeljs.io/docs/en/ if this interests you!
|
*/

const path = require("path");
const entryFile = path.resolve(__dirname, "client", "src", "index.js");
const outputDir = path.resolve(__dirname, "client", "dist");

const webpack = require("webpack");

module.exports = {
  entry: ["@babel/polyfill", entryFile],
  output: {
    path: outputDir,
    publicPath: "/",
    filename: "bundle.js",
  },
  devtool: "inline-source-map",
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        loader: "babel-loader",
        exclude: /node_modules/,
      },
      {
        test: /\.(scss|css)$/,
        use: [
          {
            loader: "style-loader",
          },
          {
            loader: "css-loader",
          },
        ],
      },
      {
        test: /\.(png|svg|jpg|gif)$/,
        use: [
          {
            loader: "url-loader",
          },
        ],
      },
    ],
  },
  resolve: {
    extensions: ["*", ".js", ".jsx"],
  },
  plugins: [new webpack.HotModuleReplacementPlugin()],
  devServer: {
    historyApiFallback: true,
    static: "./client/dist",
    hot: true,
    proxy: {
      "/api": "http://localhost:3000",
      "/socket.io/*": {
        target: "http://localhost:3000",
        ws: true,
      },
    },
  },
};
