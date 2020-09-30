

<html>
    <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Truss3D Analysis</title>
<style>
.fixedheadertable {
  width: 100%;
  max-height: 60%;
  overflow: scroll;
}

table {
  position: relative;
  border: 1px solid #ddd;
  border-collapse: collapse;
    text-decoration:none;
}

td, th {
  white-space: nowrap;
  border: 1px solid #ddd;
  padding: 5px;
  
}

th {
  background-color: #eee;
  position: -webkit-sticky;
  position: sticky;
  top: -1px;
  z-index: 2;
  text-align: center;
}
th:first-of-type {
  left: 0;
  z-index: 3;
}

tbody tr td:first-of-type {
  background-color: #eee;
  position: -webkit-sticky;
  position: sticky;
  left: -1px;
  text-align: left;
}
.hover {
  background: yellow;
}
</style>
    </head>
    <body>
<?php
include("truss3dclass.php");

$t1 = new Truss3D();

echo "<h2>Create model assigning variables</h3>";
//Create model assigning variables
$t1->Init();
echo $t1->Process();

echo "<h2>Create model from csv </h3>";
//Create model from csv 
$t1->readCSV("truss3dclass.csv");
echo $t1->Process();

?>
  
    </body>
</html>
