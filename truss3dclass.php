<?php
/******
 * Truss3DClass
 *
 * Class Inheritance
 * [FEMSolver]
 *  -> [Truss3D]
 * [FEMSolver] is base class and includes several matrix operation for the standard FEM solutions.
 * [Truss3D] is a class for the FEM solution process and include data structure of 3 dimensional pin jointed truss.
 * Solution process run for the loaded truss to analyse deformations, reactions and element forces.
 * Multiple load cases will be solved simultaneously.
 * Html result tables are generated during the process.
 * Model of 3d truss can be generated within class by assigning values to variables.
 * Or model can be created by loading CSV file.
 * Model csv file is very simple comma separated text file in-wich the properties of FEM element, boundary conditions and loads are written.
 * 
 * This class is free for the educational use as long as maintain this header together with this class.
 * Author: Win Aung Cho
 * Contact winaungcho@gmail.com
 * version 1.0
 * Date: 30-9-2020
 *
 ******/
//Starting index = 0
// relation with square matrix M[I][J] and banded Matrix BM[ib][jb]
//I=ib, J=ib+jb;
//ib=I, jb=J-ib;
class FEMSolver
{
    function bansol($s, &$f, $nq, $nbw)
    {
        /* ----- band solver ----- */
        $n1 = $nq - 1;
        /* --- forward elimination --- */
        for ($k = 1;$k <= $n1;$k++)
        {
            $nk = $nq - $k + 1;
            if ($nk > $nbw) $nk = $nbw;
            for ($i = 2;$i <= $nk;$i++)
            {
                if (!isset($s[$nbw * ($k - 1) + $i - 1])) $s[$nbw * ($k - 1) + $i - 1] = 0;
                $c1 = $s[$nbw * ($k - 1) + $i - 1] / $s[$nbw * ($k - 1) ];
                $i1 = $k + $i - 1;
                for ($j = $i;$j <= $nk;$j++)
                {
                    $j1 = $j - $i + 1;
                    if (!isset($s[$nbw * ($k - 1) + $j - 1])) $s[$nbw * ($k - 1) + $j - 1] = 0;
                    if (!isset($s[$nbw * ($i1 - 1) + $j1 - 1])) $s[$nbw * ($i1 - 1) + $j1 - 1] = 0;
                    $s[$nbw * ($i1 - 1) + $j1 - 1] = $s[$nbw * ($i1 - 1) + $j1 - 1] - $c1 * $s[$nbw * ($k - 1) + $j - 1];
                }
                if (!isset($f[$i1 - 1])) $f[$i1 - 1] = 0;
                $f[$i1 - 1] = $f[$i1 - 1] - $c1 * $f[$k - 1];
            }
        }
        /* --- back-substitution --- */
        $f[$nq - 1] = $f[$nq - 1] / $s[$nbw * ($nq - 1) ];
        for ($kk = 1;$kk <= $n1;$kk++)
        {
            $k = $nq - $kk;
            $c1 = 1.0 / $s[$nbw * ($k - 1) ];
            $f[$k - 1] = $c1 * $f[$k - 1];
            $nk = $nq - $k + 1;
            if ($nk > $nbw) $nk = $nbw;
            for ($j = 2;$j <= $nk;$j++)
            {
                $f[$k - 1] = $f[$k - 1] - $c1 * $s[$nbw * ($k - 1) + $j - 1] * $f[$k + $j - 2];
            }
        }
    }
    function Reaction($pred, $node, $cnst, $f)
    {
        $nd = 0;
        $ndn = 3;
        $n = count($pred);
        $reaction = array();
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $pred[$i]["n"];
            $reaction[$i]["nn"] = $nn;
            if ($node[$nn]["dof"][0] == 0)
            {
                $k = $ndn * $nn;
                $reaction[$i]["ux"] = $cnst * (0 - $f[$k]);
            }
            if ($node[$nn]["dof"][1] == 0)
            {
                $k = $ndn * $nn + 1;
                $reaction[$i]["uy"] = $cnst * (0 - $f[$k]);
            }
            if ($node[$nn]["dof"][2] == 0)
            {
                $k = $ndn * $nn + 2;
                $reaction[$i]["uz"] = $cnst * (0 - $f[$k]);
            }
        }
        return $reaction;
    }
    function getdispbond($pred, $node)
    {
        // u = displacement, nu = index of dof, nd = no. of boundary cond.
        $nd = 0;
        $ndn = 3;
        $n = count($pred);
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $pred[$i]["n"];
            if ($node[$nn]["dof"][0] == 0)
            {
                $u[$nd] = $pred[$i]["ux"];
                $nu[$nd] = $ndn * $nn;
                $nd++;
            }
            if ($node[$nn]["dof"][1] == 0)
            {
                $u[$nd] = $pred[$i]["uy"];
                $nu[$nd] = $ndn * $nn + 1;
                $nd++;
            }
            if ($node[$nn]["dof"][2] == 0)
            {
                $u[$nd] = $pred[$i]["uz"];
                $nu[$nd] = $ndn * $nn + 2;
                $nd++;
            }
        }
        return array(
            $nd,
            $u,
            $nu
        );
    }
    function modifydispbond($pred, $node, &$s, $nbw, $cnst, &$f)
    {
        list($nd, $u, $nu) = $this->getdispbond($pred, $node);
        // u = displacement, nu = index of dof, nd = no. of boundary cond.
        for ($i = 0;$i < $nd;$i++)
        {
            $k = $nu[$i];
            $s[$k * $nbw] = $s[$k * $nbw] + $cnst;
            if (!isset($f[$k])) $f[$k] = 0;
            $f[$k] = $f[$k] + $cnst * $u[$i];
        }
    }
    function assemble($node, $nbw, $Km, &$s)
    {
        $ndn = 3;
        $nen = 2;
        for ($ii = 0;$ii < $nen;$ii++)
        {
            $nrt = $ndn * ($node[$ii]);
            for ($it = 0;$it < $ndn;$it++)
            {
                $nr = $nrt + $it;
                $i = $ndn * $ii + $it;
                for ($jj = 0;$jj < $nen;$jj++)
                {
                    $nct = $ndn * ($node[$jj]);
                    for ($jt = 0;$jt < $ndn;$jt++)
                    {
                        $j = $ndn * $jj + $jt;
                        $nc = $nct + $jt - $nr;
                        if ($nc >= 0)
                        {
                            if (!isset($s[$nbw * $nr + $nc])) $s[$nbw * $nr + $nc] = 0;
                            $s[$nbw * $nr + $nc] += $Km[$i][$j];
                        }
                    }
                }
            }
        }
    }
    function formnodalmass(&$gM, $nbw, $nodalmass)
    //form nodal load vector
    
    {
        $ndn = 3;
        $n = count($nodalmass);
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $nodalmass[$i]["n"];

            $nnq = $nn * $ndn;
            if (!isset($gM[$nbw * $nnq])) $gM[$nbw * $nnq] = 0;
            $gM[$nbw * $nnq] = $gM[$nbw * $nnq] + $nodalmass[$i]["Mx"];

            $nnq = $nn * $ndn + 1;
            if (!isset($gM[$nbw * $nnq])) $gM[$nbw * $nnq] = 0;
            $gM[$nbw * $nnq] = $gM[$nbw * $nnq] + $nodalmass[$i]["My"];

            $nnq = $nn * $ndn + 2;
            if (!isset($gM[$nbw * $nnq])) $gM[$nbw * $nnq] = 0;
            $gM[$nbw * $nnq] = $gM[$nbw * $nnq] + $nodalmass[$i]["Mz"];
        };
    }
    function formf(&$f, $loadcase, $ln)
    //form nodal load vector
    
    {
        $ndn = 3;
        $n = count($loadcase[$ln]["nodes"]);
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $loadcase[$ln]["nodes"][$i]["n"];

            $nnq = $nn * $ndn;
            if (!isset($f[$nnq])) $f[$nnq] = 0;
            $f[$nnq] = $f[$nnq] + $loadcase[$ln]["nodes"][$i]["Fx"];

            $nnq = $nn * $ndn + 1;
            if (!isset($f[$nnq])) $f[$nnq] = 0;
            $f[$nnq] = $f[$nnq] + $loadcase[$ln]["nodes"][$i]["Fy"];

            $nnq = $nn * $ndn + 2;
            if (!isset($f[$nnq])) $f[$nnq] = 0;
            $f[$nnq] = $f[$nnq] + $loadcase[$ln]["nodes"][$i]["Fz"];
        };
    }
    function formelementload($node, $line, &$f, $loadcase, $ln)
    {
        $ndn = 3;
        if (!isset($loadcase[$ln]["dists"])) return;
        $n = count($loadcase[$ln]["dists"]);
        for ($i = 0;$i < $n;$i++)
        {
            $ne = $loadcase[$ln]["dists"][$i]["n"];
            $node1 = $node[$line[$ne]["I"]];
            $node2 = $node[$line[$ne]["J"]];
            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            $dz = $node2["z"] - $node1["z"];
            $w = array(
                $loadcase[$ln]["dists"][$i]["Wx"],
                $loadcase[$ln]["dists"][$i]["Wy"],
                $loadcase[$ln]["dists"][$i]["Wz"]
            );
            $ed = $this->elementload($w, $dx, $dy, $dz);
            $i1 = $ndn * $line[$ne]["I"];
            $i2 = $ndn * $line[$ne]["J"];

            for ($j = 0;$j < $ndn;$j++)
            {
                if (!isset($f[$i1 + $j])) $f[$i1 + $j] = 0;
                if (!isset($f[$i2 + $j])) $f[$i2 + $j] = 0;
                $f[$i1 + $j] = $f[$i1 + $j] + $ed[$j];
                $f[$i2 + $j] = $f[$i2 + $j] + $ed[$j + $ndn];
            }
        };
    }
    function KPin3D(&$Km, $EA, $dx, $dy, $dz)
    {
        // stiffness in global
        $ndn = 3;
        $nen = 2;
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2) + pow(($dz) , 2));
        $l = $dx / $ELL;
        $m = $dy / $ELL;
        $n = $dz / $ELL;
        $l2 = $l * $l;
        $m2 = $m * $m;
        $n2 = $n * $n;
        $lm = $l * $m;
        $ln = $l * $n;
        $mn = $m * $n;
        for ($i = 0;$i < $ndn * $nen;$i++) for ($j = 0;$j < $ndn * $nen;$j++)
        {
            if ($i == $j) $Km[$i][$j] = 1;
            else $Km[$i][$j] = 0.0;
        }
        $Km[0][0] = $Km[3][3] = $l2;
        $Km[1][1] = $Km[4][4] = $m2;
        $Km[2][2] = $Km[5][5] = $n2;
        $Km[0][1] = $Km[3][4] = $lm;
        $Km[0][2] = $Km[3][5] = $ln;
        $Km[1][2] = $Km[4][5] = $mn;
        $Km[0][3] = - $l2;
        $Km[1][4] = - $m2;
        $Km[2][5] = - $n2;
        $Km[0][4] = $Km[1][3] = - $lm;
        $Km[1][5] = $Km[2][4] = - $mn;
        $Km[0][5] = $Km[2][3] = - $ln;
        for ($i = 0;$i < $ndn * $nen;$i++) for ($j = $i;$j < $ndn * $nen;$j++)
        {
            $Km[$i][$j] = $Km[$i][$j] * $EA / $ELL;
            if ($i != $j) $Km[$j][$i] = $Km[$i][$j];
        }
        /*
        for ($i = 0; $i < $ndn*$nen; $i++){
        for ($j = $i; $j < $ndn*$nen; $j++)
            echo $Km[$i][$j]."::";
        echo "<br>";
        }
        */
    }
    function MPin3D(&$Mm, $mA, $dx, $dy, $dz)
    {
        // mass in global
        $ndn = 3;
        $nen = 2;
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2) + pow(($dz) , 2));
        for ($i = 0;$i < $ndn * $nen;$i++) for ($j = 0;$j < $ndn * $nen;$j++)
        {
            if ($i == $j) $Mm[$i][$j] = $mA * $ELL / 2;
            else $Mm[$i][$j] = 0.0;
        }
    }
    function elementload($w, $dx, $dy, $dz)
    {
        // element load in global, $w in global
        $ed = array();
        $ndn = 3;
        $nen = 2;
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2) + pow(($dz) , 2));
        for ($j = 0;$j < $nen;$j++)
        {
            $ed[$j * $ndn] = $w[0] * $ELL / 2;
            $ed[$j * $ndn + 1] = $w[1] * $ELL / 2;
            $ed[$j * $ndn + 2] = $w[2] * $ELL / 2;
        }
        return $ed;
    }
    function Lambada(&$L, $dx, $dy, $dz)
    {
        // transformation matrix dlocal = L dglobal
        $ndn = 3;
        $nen = 2;
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2) + pow(($dz) , 2));
        $l = $dx / $ELL;
        $m = $dy / $ELL;
        $n = $dz / $ELL;
        $L[0][0] = $l;
        $L[0][1] = $m;
        $L[0][2] = $n;
        $L[0][3] = 0;
        $L[0][4] = 0;
        $L[0][5] = 0;
        $L[1][0] = 0;
        $L[1][1] = 0;
        $L[1][2] = 0;
        $L[1][3] = $l;
        $L[1][4] = $m;
        $L[1][5] = $n;
    }
}
class Truss3D extends FEMSolver
{
    var $node = array();
    var $boundary = array();
    var $line = array();
    var $Mat = array();
    var $Sec = array();
    var $loadcase = array();
    var $nodalmass = array();
    var $processfinished = false;
    var $linecolor = "#0000ff";
    var $textcolor = "#00ff00";
    var $resulthtml = "";
    // result
    var $react;
    function Init()
    {
        $this->node[] = array(
            "x" => 0.0,
            "y" => 0.0,
            "z" => 0.0,
            "dof" => "000",
            "dx" => 0.00,
            "dy" => 0.00,
            "dz" => 0.0
        );
        $this->node[] = array(
            "x" => 40.0,
            "y" => 0.0,
            "z" => 0.0,
            "dof" => "100",
            "dx" => 0.0,
            "dy" => 0.0,
            "dz" => 0.0
        );
        $this->node[] = array(
            "x" => 40.0,
            "y" => 30.0,
            "z" => 0.0,
            "dof" => "110",
            "dx" => 0.00,
            "dy" => 0.0,
            "dz" => 0.0
        );
        $this->node[] = array(
            "x" => 0.0,
            "y" => 30.0,
            "z" => 0.0,
            "dof" => "000",
            "dx" => 0.00,
            "dy" => 0.00,
            "dz" => 0.0
        );
        $this->line[] = array(
            "I" => 0,
            "J" => 1,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->line[] = array(
            "I" => 1,
            "J" => 2,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->line[] = array(
            "I" => 2,
            "J" => 3,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->line[] = array(
            "I" => 0,
            "J" => 2,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->boundary[] = array(
            "n" => 0,
            "ux" => 0.0,
            "uy" => 0.0,
            "uz" => 0.0
        );
        $this->boundary[] = array(
            "n" => 1,
            "uy" => 0.0,
            "uz" => 0.0
        );
        $this->boundary[] = array(
            "n" => 2,
            "uz" => 0.0
        );
        $this->boundary[] = array(
            "n" => 3,
            "ux" => 0.0,
            "uy" => 0.0,
            "uz" => 0.0
        );
        $this->Mat[] = array(
            "name" => "Steel",
            "E" => 29.500e6,
            "alpha" => 6.5e6,
            "wpv" => 0.2836,
            "mpv" => 7.345E-04
        );
        $this->Sec[] = array(
            "name" => "L2x2",
            "A" => 1.0,
            "I" => 10
        );
        $this->Sec[] = array(
            "name" => "L3x3",
            "A" => 2.0,
            "I" => 20
        );
        $this->nodalmass[] = array(
            "n" => 1,
            "Mx" => 0.0,
            "My" => 0.0,
            "Mz" => 0.0
        );
        $this->loadcase[0]["nodes"][] = array(
            "n" => 1,
            "Fx" => 20000.0,
            "Fy" => 0.0,
            "Fz" => 0.0
        );
        $this->loadcase[0]["nodes"][] = array(
            "n" => 2,
            "Fx" => 0.0,
            "Fy" => - 25000.0,
            "Fz" => 0.0
        );
    }
    function Reset()
    {
        $this->node = array();
        $this->line = array();
        $this->boundary = array();
        $this->Mat = array();
        $this->Sec = array();
        $this->nodalmass = array();
        $this->loadcase = array();
    }
    function Process() // analyse on frame
    
    {
        $ndn = 3;
        $nen = 2;
        $itmax = 50;
        $s = array();
        $f = array();
        $gM = array();
        if (!$this->node) return "No valid Pin3D frame data found";
        $nn = count($this->node);
        $nq = $nn * $ndn;
        $ne = count($this->line);
        $nbw = 6;

        for ($i = 0;$i < $ne;$i++)
        {
            $nb = $ndn * (abs($this->line[$i]["I"] - $this->line[$i]["J"]) + 1.0);
            if ($nbw < $nb) $nbw = $nb;
        }
        for ($i = 0;$i < $nq;$i++)
        {
            $st[$i] = rand(-10, 10) / 10.0;
            for ($j = 0;$j < $nbw;$j++)
            {
                $s[$nbw * $i + $j] = 0;
                $gM[$nbw * $i + $j] = 0;
            }
        }
        //print("NBW=".$nbw);
        for ($i = 0;$i < $ne;$i++)
        {
            $nmat = $this->line[$i]["Mat"];
            $nsec = $this->line[$i]["Sec"];
            $node1 = $this->node[$this->line[$i]["I"]];
            $node2 = $this->node[$this->line[$i]["J"]];
            $EA = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["A"];
            $mA = $this->Mat[$nmat]["mpv"] * $this->Sec[$nsec]["A"];
            //echo "Mass:".$mA."<br>";
            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            $dz = $node2["z"] - $node1["z"];
            $Km = array();
            $Mm = array();
            $this->KPin3D($Km, $EA, $dx, $dy, $dz);
            $node = array(
                $this->line[$i]["I"],
                $this->line[$i]["J"]
            );
            $this->assemble($node, $nbw, $Km, $s);
            $this->MPin3D($Mm, $mA, $dx, $dy, $dz);
            $this->assemble($node, $nbw, $Mm, $gM);
            //printf("Element %d <br>", $i);
            
        }

        $this->formnodalmass($gM, $nbw, $this->nodalmass);
        $cnst = 0;
        for ($i = 0;$i < $nq;$i++)
        {
            if ($cnst < $s[$i * $nbw]) $cnst = $s[$i * $nbw];
        }
        $cnst = 10000 * $cnst;

        $n = count($this->loadcase);
        for ($i = 0;$i < $nq;$i++) $f[$i] = 0;
        for ($i = 0;$i < $n;$i++)
        {
            $this->formf($f, $this->loadcase, $i);
            $this->formelementload($this->node, $this->line, $f, $this->loadcase, $i);
        }
        $this->modifydispbond($this->boundary, $this->node, $s, $nbw, $cnst, $f);
        $s2 = $s;
        $sh = 0;
        $nev = 3;

        $this->bansol($s, $f, $nq, $nbw);

        // displacement table
        $outdisp = "<h4>Nodal Displacements</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\"><tr><th>node#</th><th>x-displ.</th><th>y-displ.</th><th>z-displ</th></tr>";
        $nn = count($this->node);
        for ($i = 0;$i < $nn;$i++)
        {
            $i1 = $ndn * $i;
            $str = sprintf("<tr><td align='right'>%3d</td><td align='right'>%01.8f</td><td align='right'>%01.8f</td><td align='right'>%01.8f</td></tr>", $i, $f[$i1], $f[$i1 + 1], $f[$i1 + 2]);
            $outdisp .= $str;
            $this->node[$i]["dx"] = $f[$i1];
            $this->node[$i]["dy"] = $f[$i1 + 1];
            $this->node[$i]["dz"] = $f[$i1 + 2];
        }
        $outdisp .= "</table></div>";

        $ne = count($this->line);
        $outfrc = "<h3>Element Fources</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\"><tr><th>member#</th><th>XI</th><th>YI</th><th>ZI</th><th>XJ</th><th>YJ</th><th>ZJ</th></tr>";
        for ($i = 0;$i < $ne;$i++)
        {
            $nmat = $this->line[$i]["Mat"];
            $nsec = $this->line[$i]["Sec"];
            $node1 = $this->node[$this->line[$i]["I"]];
            $node2 = $this->node[$this->line[$i]["J"]];
            $EA = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["A"];

            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            $dz = $node2["z"] - $node1["z"];
            $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2) + pow(($dz) , 2));
            //MessageOut(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $EA, $EI, $dx, $dy);
            //MessageOut(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $this->Mat[$nmat]["E"], $this->Sec[$nsec]["I"], $dx, $dy);
            $Km = array();
            $this->KPin3D($Km, $EA, $dx, $dy, $dz);
            list($i1, $i2) = array(
                $this->line[$i]["I"] * $ndn,
                $this->line[$i]["J"] * $ndn
            );
            //$alambada = lambada($dx, $dy, $dz);
            $ed = array();
            for ($j = 0;$j < $ndn;$j++)
            {
                $ed[$j] = $f[$i1 + $j];
                $ed[$j + $ndn] = $f[$i2 + $j];
            }
            $edp = array();
            // dp = Ld
            for ($j = 0;$j < $nen * $ndn;$j++)
            {
                //$edp[$j] = 0;
                //for ($k = 0; $k < $nen*$ndn; $k++) {
                //    $edp[$j] = $edp[$j] + $alambada[$j][$k] * $ed[$k];
                //}
                $edp[$j] = $ed[$j];
            }
            $ed = array();
            for ($j = 0;$j < $nen * $ndn;$j++) $ed[$j] = 0;
            // $ed should member load vector
            $nl = count($this->loadcase);

            for ($j = 0;$j < $nl;$j++)
            {
                if (!isset($this->loadcase[$j]['dists'])) continue;
                $nk = count($this->loadcase[$j]['dists']);
                for ($k = 0;$k < $nk;$k++)
                {
                    $mn = $this->loadcase[$j]['dists'][$k]["n"];
                    if ($i == $mn)
                    {
                        $w = array(
                            $this->loadcase[$j]["dists"][$k]["Wx"],
                            $this->loadcase[$j]["dists"][$k]["Wy"],
                            $this->loadcase[$j]["dists"][$k]["Wz"]
                        );
                        $ed = $this->elementload($w, $dx, $dy, $dz);
                    }
                }
            }
            for ($j = 0;$j < $nen * $ndn;$j++)
            {
                $ef[$j] = - $ed[$j];
                for ($k = 0;$k < $nen * $ndn;$k++)
                {
                    $ef[$j] = $ef[$j] + $Km[$j][$k] * $edp[$k]; // $Km must be global
                    
                }
            }
            $str = sprintf("<tr ><td>%d", $i);
            $outfrc .= $str;
            for ($j = 0;$j < 2;$j++)
            {
                $ii = $ndn * $j;
                $str = sprintf("</td><td align='right'>%11.4f</td><td align='right'>%11.4f</td><td align='right'>%11.4f", $ef[$ii], $ef[$ii + 1], $ef[$ii + 2]);
                $outfrc .= $str;
            }
            $outfrc .= "</td></tr>";
        }
        $outfrc .= "</table></div>";

        $this->react = $this->Reaction($this->boundary, $this->node, $cnst, $f);
        $n = count($this->react);
        $outreac = "<h3>Reaction</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\"><tr><th>node#</th><th>Rx</th><th>Ry</th><th>Rz</th></tr>";
        for ($i = 0;$i < $n;$i++)
        {
            $rx = $ry = $dz = 0;
            if (isset($this->react[$i]["ux"])) $rx = $this->react[$i]["ux"];
            if (isset($this->react[$i]["uy"])) $ry = $this->react[$i]["uy"];
            if (isset($this->react[$i]["uz"])) $rz = $this->react[$i]["uz"];
            $str = sprintf("<tr><td align='right'>%3d</td><td align='right'>%01.4f</td><td align='right'>%01.4f</td><td align='right'>%01.4f</td></tr>", $this->react[$i]["nn"], $rx, $ry, $rz);
            $outreac .= $str;
        };
        $outreac .= "</table>";
        $outreac .= "</div>";

        $this->resulthtml = $outdisp . $outfrc . $outreac;
        $this->processfinished = true;
        return $this->resulthtml;
    }
    function readCSV($fname)
    {
        $this->Reset();
        $handle = fopen($fname, "r");
        if ($handle)
        {
            while (!feof($handle) && ($line = fgets($handle)) !== false)
            {
                // process the line read.
                if (!empty($line) && $line[0] !== ';')
                {
                    $args = explode(",", $line);
                    switch ($args[0])
                    {
                        case "node":
                            $this->node[] = array(
                                "x" => (float)$args[1],
                                "y" => (float)$args[2],
                                "z" => (float)$args[3],
                                "dof" => $args[4],
                                "dx" => 0.00,
                                "dy" => 0.00,
                                "dz" => 0.0
                            );
                        break;
                        case "line":
                            $this->line[] = array(
                                "I" => (int)$args[1],
                                "J" => (int)$args[2],
                                "Mat" => (int)$args[3],
                                "Sec" => (int)$args[4]
                            );
                        break;
                        case "boundary":
                            $this->boundary[] = array(
                                "n" => (int)$args[1],
                                "ux" => (float)$args[2],
                                "uy" => (float)$args[3],
                                "uz" => (float)$args[4]
                            );
                        break;
                        case "Mat":
                            $this->Mat[] = array(
                                "name" => $args[1],
                                "E" => (float)$args[2],
                                "alpha" => (float)$args[3],
                                "wpv" => (float)$args[4],
                                "mpv" => (float)$args[5]
                            );
                        break;
                        case "Sec":
                            $this->Sec[] = array(
                                "name" => $args[1],
                                "A" => (float)$args[2],
                                "I" => (float)$args[3]
                            );
                        break;
                        case "nodalmass":
                            $this->nodalmass[] = array(
                                "n" => (int)$args[1],
                                "Mx" => (float)$args[2],
                                "My" => (float)$args[3],
                                "Mz" => (float)$args[4]
                            );
                        break;
                        case "loadcase":
                            if ($args[2] === "nodes") $this->loadcase[(int)$args[1]]["nodes"][] = array(
                                "n" => (int)$args[3],
                                "Fx" => (float)$args[4],
                                "Fy" => (float)$args[5],
                                "Fz" => (float)$args[6]
                            );
                            if ($args[2] === "dists") $this->loadcase[(int)$args[1]]["dists"][] = array(
                                "n" => (int)$args[3],
                                "Wx" => (float)$args[4],
                                "Wy" => (float)$args[5],
                                "Wz" => (float)$args[6]
                            );
                            break;
                        }
                    }

            }
            fclose($handle);
        }
        else
        {
            // error opening the file.
            
        }
    }
}

?>
