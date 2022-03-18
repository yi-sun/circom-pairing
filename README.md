# circom-pairing

TODO: add table of contents with links

## Setup
First, install [yarn](https://classic.yarnpkg.com/en/) and [circom](https://docs.circom.io/getting-started/installation/). 
Then run ```yarn install``` in the root directory to install the dependencies in ```yarn.lock```. 

## Testing
See the ```/test``` directory for examples of tests. The circuits to be tested should be written in the ```/test/circuits``` folder, while the test execution code should be written in regular JavaScript files under ```/test```. A short description of each test can be passed in as the first parameter of the ```describe()``` function, and ```yarn --grep name``` will run all tests whose description contains ```name``` as a substring. 