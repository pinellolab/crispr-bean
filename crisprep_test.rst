Tutorial
==================

Load the required packages. (Anndata import isn't required to use the package).

.. code:: ipython3

    import numpy as np
    import pandas as pd
    import anndata as ad
    import seaborn as sns
    import matplotlib.pyplot as plt
    import beret as cp


beret `ReporterScreen` object and perturb-seq Screen` object are both `anndata` compatible.

.. code:: ipython3

    adata = ad.read_h5ad("beret_count_07+1021_LDLvar.h5ad")

.. code:: ipython3

    adata


.. parsed-literal::

    AnnData object with n_obs × n_vars = 3455 × 12
        obs: 'name', 'Unnamed: 0', 'Target gene/variant', 'Target descriptor', 'Arbitrary number', 'gRNA position category', 'Target base position in gRNA', 'Target base position in reporter', 'BE', 'Group', 'sequence', 'Reporter', 'barcode', '5-nt PAM', 'offset', 'target', 'target_pos', 'Group2', 'masked_sequence', 'masked_barcode', 'edit_rate'
        var: 'index', 'sort', 'replicate'
        uns: 'allele_counts', 'edit_counts'
        layers: 'X_bcmatch', 'edits'



.. code:: ipython3

    cdata = cp.read_h5ad("beret_count_07+1021_LDLvar.h5ad")

.. code:: ipython3

    cdata





.. parsed-literal::

    Genome Editing Screen comprised of n_guides x n_conditions = 3455 x 12
       guides:    'name', 'Unnamed: 0', 'Target gene/variant', 'Target descriptor', 'Arbitrary number', 'gRNA position category', 'Target base position in gRNA', 'Target base position in reporter', 'BE', 'Group', 'sequence', 'Reporter', 'barcode', '5-nt PAM', 'offset', 'target', 'target_pos', 'Group2', 'masked_sequence', 'masked_barcode', 'edit_rate'
       condit:    'index', 'sort', 'replicate'
       condit_m:  
       condit_p:  
       layers:    'X_bcmatch', 'edits'
       uns:       'allele_counts', 'edit_counts'

-  cdata.X: guide count
-  cdata.guides: guide metadata
-  cdata.condit: sample/condition metadata
-  cdata.layers["X_bcmatch"]: barcode-matched guide counts
-  cdata.layers["edits"]: edit counts
-  cdata.uns["allele_counts"]: allele counts per guide and condition
-  cdata.uns["edit_counts"]: edit counts per guide and condition

`guides` attribute contains the information about each guide.

.. code:: ipython3

    cdata.guides





.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>name</th>
          <th>Unnamed: 0</th>
          <th>Target gene/variant</th>
          <th>Target descriptor</th>
          <th>Arbitrary number</th>
          <th>gRNA position category</th>
          <th>Target base position in gRNA</th>
          <th>Target base position in reporter</th>
          <th>BE</th>
          <th>Group</th>
          <th>...</th>
          <th>Reporter</th>
          <th>barcode</th>
          <th>5-nt PAM</th>
          <th>offset</th>
          <th>target</th>
          <th>target_pos</th>
          <th>Group2</th>
          <th>masked_sequence</th>
          <th>masked_barcode</th>
          <th>edit_rate</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CONTROL_1_g1</td>
          <td>0</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g1</td>
          <td>4</td>
          <td>10</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>CCAAGCCCTACGCGGTAGGGAACTTTGGGAGC</td>
          <td>GTTT</td>
          <td>GGGAG</td>
          <td>-10</td>
          <td>CONTROL_1</td>
          <td>9</td>
          <td>NegCtrl</td>
          <td>CCTGCGCGGTGGGGGGCTTT</td>
          <td>GTTT</td>
          <td>0.531163</td>
        </tr>
        <tr>
          <th>1</th>
          <td>CONTROL_1_g2</td>
          <td>1</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g2</td>
          <td>5</td>
          <td>11</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>TCCAAGCCCTACGCGGTAGGGAACTTTGGGAG</td>
          <td>AACA</td>
          <td>TGGGA</td>
          <td>-11</td>
          <td>CONTROL_1</td>
          <td>10</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCTT</td>
          <td>GGCG</td>
          <td>0.640765</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CONTROL_1_g3</td>
          <td>2</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g3</td>
          <td>5</td>
          <td>12</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>GTCCAAGCCCTACGCGGTAGGGAACTTTGGGA</td>
          <td>CGCT</td>
          <td>TTGGG</td>
          <td>-12</td>
          <td>CONTROL_1</td>
          <td>11</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCT</td>
          <td>CGCT</td>
          <td>0.417709</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CONTROL_1_g4</td>
          <td>3</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>CGTCCAAGCCCTACGCGGTAGGGAACTTTGGG</td>
          <td>TGAG</td>
          <td>TTTGG</td>
          <td>-13</td>
          <td>CONTROL_1</td>
          <td>12</td>
          <td>NegCtrl</td>
          <td>GGCCCTGCGCGGTGGGGGGC</td>
          <td>TGGG</td>
          <td>0.126400</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CONTROL_1_g5</td>
          <td>4</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>ACGTCCAAGCCCTACGCGGTAGGGAACTTTGG</td>
          <td>GTAT</td>
          <td>CTTTG</td>
          <td>-14</td>
          <td>CONTROL_1</td>
          <td>13</td>
          <td>NegCtrl</td>
          <td>GGGCCCTGCGCGGTGGGGGG</td>
          <td>GTGT</td>
          <td>0.201104</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>3450</th>
          <td>rs9987289_Maj_ABE_347_g1</td>
          <td>3450</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g1</td>
          <td>3</td>
          <td>10</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TGCTTGGGCATCAATATCACGTGGAACCAGCC</td>
          <td>CAGT</td>
          <td>CCAGC</td>
          <td>-10</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>9</td>
          <td>Variant</td>
          <td>GCGTCGGTGTCGCGTGGGG</td>
          <td>CGGT</td>
          <td>0.087379</td>
        </tr>
        <tr>
          <th>3451</th>
          <td>rs9987289_Maj_ABE_347_g2</td>
          <td>3451</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g2</td>
          <td>4</td>
          <td>11</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>ATGCTTGGGCATCAATATCACGTGGAACCAGC</td>
          <td>TCGC</td>
          <td>ACCAG</td>
          <td>-11</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>10</td>
          <td>Variant</td>
          <td>GGCGTCGGTGTCGCGTGGG</td>
          <td>TCGC</td>
          <td>0.299923</td>
        </tr>
        <tr>
          <th>3452</th>
          <td>rs9987289_Maj_ABE_347_g3</td>
          <td>3452</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g3</td>
          <td>6</td>
          <td>12</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>GATGCTTGGGCATCAATATCACGTGGAACCAG</td>
          <td>GCAC</td>
          <td>AACCA</td>
          <td>-12</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>11</td>
          <td>Variant</td>
          <td>TGGGCGTCGGTGTCGCGTGG</td>
          <td>GCGC</td>
          <td>0.224973</td>
        </tr>
        <tr>
          <th>3453</th>
          <td>rs9987289_Maj_ABE_347_g4</td>
          <td>3453</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>AGATGCTTGGGCATCAATATCACGTGGAACCA</td>
          <td>TTGC</td>
          <td>GAACC</td>
          <td>-13</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>12</td>
          <td>Variant</td>
          <td>TTGGGCGTCGGTGTCGCGTG</td>
          <td>TTGC</td>
          <td>0.265378</td>
        </tr>
        <tr>
          <th>3454</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>3454</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TAGATGCTTGGGCATCAATATCACGTGGAACC</td>
          <td>GCGA</td>
          <td>GGAAC</td>
          <td>-14</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>13</td>
          <td>Variant</td>
          <td>CTTGGGCGTCGGTGTCGCGT</td>
          <td>GCGG</td>
          <td>0.266573</td>
        </tr>
      </tbody>
    </table>
    <p>3455 rows × 21 columns</p>
    </div>


`condit` attribute contains the sample and condition specific information.

.. code:: ipython3

    cdata.condit





.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>index</th>
          <th>sort</th>
          <th>replicate</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>rep1_bot</td>
          <td>bot</td>
          <td>rep1</td>
        </tr>
        <tr>
          <th>1</th>
          <td>rep2_bot</td>
          <td>bot</td>
          <td>rep2</td>
        </tr>
        <tr>
          <th>2</th>
          <td>rep3_VPA_bot</td>
          <td>bot</td>
          <td>rep3_VPA</td>
        </tr>
        <tr>
          <th>3</th>
          <td>rep4_VPA_bot</td>
          <td>bot</td>
          <td>rep4_VPA</td>
        </tr>
        <tr>
          <th>4</th>
          <td>rep1_bulk</td>
          <td>bulk</td>
          <td>rep1</td>
        </tr>
        <tr>
          <th>5</th>
          <td>rep2_bulk</td>
          <td>bulk</td>
          <td>rep2</td>
        </tr>
        <tr>
          <th>6</th>
          <td>rep3_VPA_bulk</td>
          <td>bulk</td>
          <td>rep3_VPA</td>
        </tr>
        <tr>
          <th>7</th>
          <td>rep4_VPA_bulk</td>
          <td>bulk</td>
          <td>rep4_VPA</td>
        </tr>
        <tr>
          <th>8</th>
          <td>rep1_top</td>
          <td>top</td>
          <td>rep1</td>
        </tr>
        <tr>
          <th>9</th>
          <td>rep2_top</td>
          <td>top</td>
          <td>rep2</td>
        </tr>
        <tr>
          <th>10</th>
          <td>rep3_VPA_top</td>
          <td>top</td>
          <td>rep3_VPA</td>
        </tr>
        <tr>
          <th>11</th>
          <td>rep4_VPA_top</td>
          <td>top</td>
          <td>rep4_VPA</td>
        </tr>
      </tbody>
    </table>
    </div>


Allele_counts information is stored in `.uns["allele_counts"]`.

.. code:: ipython3

    cdata.uns["allele_counts"]





.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>guide</th>
          <th>allele</th>
          <th>rep1_bot</th>
          <th>rep2_bot</th>
          <th>rep3_VPA_bot</th>
          <th>rep4_VPA_bot</th>
          <th>rep1_bulk</th>
          <th>rep2_bulk</th>
          <th>rep3_VPA_bulk</th>
          <th>rep4_VPA_bulk</th>
          <th>rep1_top</th>
          <th>rep2_top</th>
          <th>rep3_VPA_top</th>
          <th>rep4_VPA_top</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>0:9:+:A&gt;G,5:14:+:A&gt;G</td>
          <td>14</td>
          <td>20</td>
          <td>13</td>
          <td>0</td>
          <td>6</td>
          <td>15</td>
          <td>2</td>
          <td>17</td>
          <td>22</td>
          <td>14</td>
          <td>34</td>
          <td>3</td>
        </tr>
        <tr>
          <th>1</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-4:5:+:A&gt;G,-2:7:+:A&gt;G,5:14:+:A&gt;G,10:19:+:A&gt;G</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>2</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-7:2:+:A&gt;G,0:9:+:A&gt;G,5:14:+:A&gt;G</td>
          <td>3</td>
          <td>4</td>
          <td>2</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>5</td>
          <td>2</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
        </tr>
        <tr>
          <th>3</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-9:0:+:G&gt;A,-8:1:+:G&gt;A,-7:2:+:A&gt;C,-6:3:+:C&gt;A,-4...</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>2</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
        </tr>
        <tr>
          <th>4</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-7:2:+:A&gt;G,10:19:+:A&gt;G</td>
          <td>1</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>438407</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>4:17:+:A&gt;G,6:19:+:A&gt;G,9:22:+:A&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>2</td>
          <td>0</td>
        </tr>
        <tr>
          <th>438408</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>-12:1:+:A&gt;G,6:19:+:A&gt;G,9:22:+:A&gt;G,11:24:+:G&gt;A</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
        </tr>
        <tr>
          <th>438409</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>-12:1:+:A&gt;G,6:19:+:A&gt;G,9:22:+:A&gt;G,16:29:+:A&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
        </tr>
        <tr>
          <th>438410</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>-12:1:+:A&gt;G,0:13:+:A&gt;G,6:19:+:A&gt;G,9:22:+:A&gt;G,1...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>438411</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>-12:1:+:A&gt;G,6:19:+:A&gt;G,9:22:+:A&gt;G,12:25:+:T&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
        </tr>
      </tbody>
    </table>
    <p>438412 rows × 14 columns</p>
    </div>


Base-level edit counts can be saved at `.uns["edit_counts"]`.

.. code:: ipython3

    cdata.uns["edit_counts"]





.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>guide</th>
          <th>edit</th>
          <th>rep1_bot</th>
          <th>rep2_bot</th>
          <th>rep3_VPA_bot</th>
          <th>rep4_VPA_bot</th>
          <th>rep1_bulk</th>
          <th>rep2_bulk</th>
          <th>rep3_VPA_bulk</th>
          <th>rep4_VPA_bulk</th>
          <th>rep1_top</th>
          <th>rep2_top</th>
          <th>rep3_VPA_top</th>
          <th>rep4_VPA_top</th>
          <th>ref_base</th>
          <th>alt_base</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-1:8:+:G&gt;A</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>G</td>
          <td>A</td>
        </tr>
        <tr>
          <th>1</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-1:8:+:G&gt;C</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>G</td>
          <td>C</td>
        </tr>
        <tr>
          <th>2</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-1:8:+:G&gt;T</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>G</td>
          <td>T</td>
        </tr>
        <tr>
          <th>3</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-2:7:+:A&gt;C</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>2</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>A</td>
          <td>C</td>
        </tr>
        <tr>
          <th>4</th>
          <td>12:51779544AGA_Maj_ABE_2_g1</td>
          <td>-2:7:+:A&gt;G</td>
          <td>19</td>
          <td>34</td>
          <td>40</td>
          <td>4</td>
          <td>59</td>
          <td>25</td>
          <td>66</td>
          <td>7</td>
          <td>68</td>
          <td>48</td>
          <td>149</td>
          <td>2</td>
          <td>A</td>
          <td>G</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>217563</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>8:21:+:C&gt;A</td>
          <td>0</td>
          <td>7</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>1</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>C</td>
          <td>A</td>
        </tr>
        <tr>
          <th>217564</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>8:21:+:C&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>2</td>
          <td>0</td>
          <td>0</td>
          <td>8</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>1</td>
          <td>8</td>
          <td>0</td>
          <td>C</td>
          <td>G</td>
        </tr>
        <tr>
          <th>217565</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>8:21:+:C&gt;T</td>
          <td>0</td>
          <td>0</td>
          <td>7</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>7</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>C</td>
          <td>T</td>
        </tr>
        <tr>
          <th>217566</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>9:22:+:A&gt;G</td>
          <td>9</td>
          <td>21</td>
          <td>30</td>
          <td>51</td>
          <td>37</td>
          <td>46</td>
          <td>12</td>
          <td>20</td>
          <td>58</td>
          <td>23</td>
          <td>59</td>
          <td>47</td>
          <td>A</td>
          <td>G</td>
        </tr>
        <tr>
          <th>217567</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>9:22:+:A&gt;T</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>7</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>A</td>
          <td>T</td>
        </tr>
      </tbody>
    </table>
    <p>217568 rows × 16 columns</p>
    </div>





Subsetting & addition
---------------------

Works as anndata, supports allele & edit count operations.

Subsetting & selection
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    cdata_subset = cdata[:10,cdata.condit.sort == "bulk"]


.. parsed-literal::

    ['rep1_bulk', 'rep2_bulk', 'rep3_VPA_bulk', 'rep4_VPA_bulk']


.. code:: ipython3

    cdata_subset.uns["allele_counts"]




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>guide</th>
          <th>allele</th>
          <th>rep1_bulk</th>
          <th>rep2_bulk</th>
          <th>rep3_VPA_bulk</th>
          <th>rep4_VPA_bulk</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>14979</th>
          <td>CONTROL_10_g1</td>
          <td>-4:5:+:A&gt;G,0:9:+:A&gt;G</td>
          <td>8</td>
          <td>1</td>
          <td>3</td>
          <td>0</td>
        </tr>
        <tr>
          <th>14980</th>
          <td>CONTROL_10_g1</td>
          <td>-7:2:+:C&gt;T</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>10</td>
        </tr>
        <tr>
          <th>14981</th>
          <td>CONTROL_10_g1</td>
          <td>-4:5:+:A&gt;G</td>
          <td>29</td>
          <td>2</td>
          <td>29</td>
          <td>25</td>
        </tr>
        <tr>
          <th>14982</th>
          <td>CONTROL_10_g1</td>
          <td>1:10:+:A&gt;G</td>
          <td>0</td>
          <td>6</td>
          <td>4</td>
          <td>1</td>
        </tr>
        <tr>
          <th>14983</th>
          <td>CONTROL_10_g1</td>
          <td>-4:5:+:A&gt;G,1:10:+:A&gt;G</td>
          <td>1</td>
          <td>11</td>
          <td>5</td>
          <td>12</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>22837</th>
          <td>CONTROL_1_g5</td>
          <td>-13:0:+:A&gt;-,-12:1:+:C&gt;T,-9:4:+:C&gt;G,-8:5:+:C&gt;T,...</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>22838</th>
          <td>CONTROL_1_g5</td>
          <td>-6:7:+:A&gt;C,7:20:+:A&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>22839</th>
          <td>CONTROL_1_g5</td>
          <td>-13:0:+:A&gt;G,-10:3:+:T&gt;G,0:13:+:A&gt;G,7:20:+:A&gt;G</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>22840</th>
          <td>CONTROL_1_g5</td>
          <td>0:13:+:A&gt;T</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>22841</th>
          <td>CONTROL_1_g5</td>
          <td>0:13:+:A&gt;G,18:31:+:G&gt;A</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    <p>1080 rows × 6 columns</p>
    </div>



LFC calculation & Addition
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    cdata1 = cp.read_h5ad("/data/pinello/PROJECTS/2021_08_ANBE/data/072121_ABE_topbot/beret_counts/LDLvar/032422_crispresso/beret_count_072121_ABE_topbot_LDLvar.h5ad")
    cdata2 = cp.read_h5ad("/data/pinello/PROJECTS/2021_08_ANBE/data/102121_ABE_topbot/beret_counts/LDLvar/032422_crispresso/beret_count_102121_ABE_topbot_LDLvar.h5ad")


.. code:: ipython3

    cdata1.condit["sort"] = cdata1.condit["index"].map(lambda s: s.rsplit("_", 1)[-1])
    cdata1.condit["replicate"] = cdata1.condit["index"].map(lambda s: s.rsplit("_", 1)[0])
    cdata2.condit["sort"] = cdata2.condit["index"].map(lambda s: s.rsplit("_", 1)[-1])
    cdata2.condit["replicate"] = cdata2.condit["index"].map(lambda s: s.rsplit("_", 1)[0])

.. code:: ipython3

    cdata1.log_norm()
    lfc1 = cdata1.log_fold_change_reps("bot", "top")
    cdata2.log_norm()
    lfc2 = cdata2.log_fold_change_reps("bot", "top")
    lfcs = lfc1.join(lfc2, lsuffix = "_1", rsuffix = "_2")
    sns.pairplot(lfcs)


.. image:: output_20_2.png


LFC can be aggregated for biological replicates.

.. code:: ipython3

    cdata1.log_fold_change_aggregate("bot", "top", aggregate_condit = "replicate")

.. code:: ipython3

    cdata1.guides




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>name</th>
          <th>Unnamed: 0</th>
          <th>Target gene/variant</th>
          <th>Target descriptor</th>
          <th>Arbitrary number</th>
          <th>gRNA position category</th>
          <th>Target base position in gRNA</th>
          <th>Target base position in reporter</th>
          <th>BE</th>
          <th>Group</th>
          <th>...</th>
          <th>Reporter</th>
          <th>barcode</th>
          <th>5-nt PAM</th>
          <th>offset</th>
          <th>target</th>
          <th>target_pos</th>
          <th>Group2</th>
          <th>masked_sequence</th>
          <th>masked_barcode</th>
          <th>bot_top.lfc.median</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CONTROL_1_g1</td>
          <td>0</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g1</td>
          <td>4</td>
          <td>10</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>CCAAGCCCTACGCGGTAGGGAACTTTGGGAGC</td>
          <td>GTTT</td>
          <td>GGGAG</td>
          <td>-10</td>
          <td>CONTROL_1</td>
          <td>9</td>
          <td>NegCtrl</td>
          <td>CCTGCGCGGTGGGGGGCTTT</td>
          <td>GTTT</td>
          <td>-0.158787</td>
        </tr>
        <tr>
          <th>1</th>
          <td>CONTROL_1_g2</td>
          <td>1</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g2</td>
          <td>5</td>
          <td>11</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>TCCAAGCCCTACGCGGTAGGGAACTTTGGGAG</td>
          <td>AACA</td>
          <td>TGGGA</td>
          <td>-11</td>
          <td>CONTROL_1</td>
          <td>10</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCTT</td>
          <td>GGCG</td>
          <td>-0.212254</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CONTROL_1_g3</td>
          <td>2</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g3</td>
          <td>5</td>
          <td>12</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>GTCCAAGCCCTACGCGGTAGGGAACTTTGGGA</td>
          <td>CGCT</td>
          <td>TTGGG</td>
          <td>-12</td>
          <td>CONTROL_1</td>
          <td>11</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCT</td>
          <td>CGCT</td>
          <td>0.186679</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CONTROL_1_g4</td>
          <td>3</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>CGTCCAAGCCCTACGCGGTAGGGAACTTTGGG</td>
          <td>TGAG</td>
          <td>TTTGG</td>
          <td>-13</td>
          <td>CONTROL_1</td>
          <td>12</td>
          <td>NegCtrl</td>
          <td>GGCCCTGCGCGGTGGGGGGC</td>
          <td>TGGG</td>
          <td>-0.022441</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CONTROL_1_g5</td>
          <td>4</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>ACGTCCAAGCCCTACGCGGTAGGGAACTTTGG</td>
          <td>GTAT</td>
          <td>CTTTG</td>
          <td>-14</td>
          <td>CONTROL_1</td>
          <td>13</td>
          <td>NegCtrl</td>
          <td>GGGCCCTGCGCGGTGGGGGG</td>
          <td>GTGT</td>
          <td>0.457033</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>3450</th>
          <td>rs9987289_Maj_ABE_347_g1</td>
          <td>3450</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g1</td>
          <td>3</td>
          <td>10</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TGCTTGGGCATCAATATCACGTGGAACCAGCC</td>
          <td>CAGT</td>
          <td>CCAGC</td>
          <td>-10</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>9</td>
          <td>Variant</td>
          <td>GCGTCGGTGTCGCGTGGGG</td>
          <td>CGGT</td>
          <td>-0.418312</td>
        </tr>
        <tr>
          <th>3451</th>
          <td>rs9987289_Maj_ABE_347_g2</td>
          <td>3451</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g2</td>
          <td>4</td>
          <td>11</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>ATGCTTGGGCATCAATATCACGTGGAACCAGC</td>
          <td>TCGC</td>
          <td>ACCAG</td>
          <td>-11</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>10</td>
          <td>Variant</td>
          <td>GGCGTCGGTGTCGCGTGGG</td>
          <td>TCGC</td>
          <td>-0.084936</td>
        </tr>
        <tr>
          <th>3452</th>
          <td>rs9987289_Maj_ABE_347_g3</td>
          <td>3452</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g3</td>
          <td>6</td>
          <td>12</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>GATGCTTGGGCATCAATATCACGTGGAACCAG</td>
          <td>GCAC</td>
          <td>AACCA</td>
          <td>-12</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>11</td>
          <td>Variant</td>
          <td>TGGGCGTCGGTGTCGCGTGG</td>
          <td>GCGC</td>
          <td>-0.339419</td>
        </tr>
        <tr>
          <th>3453</th>
          <td>rs9987289_Maj_ABE_347_g4</td>
          <td>3453</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>AGATGCTTGGGCATCAATATCACGTGGAACCA</td>
          <td>TTGC</td>
          <td>GAACC</td>
          <td>-13</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>12</td>
          <td>Variant</td>
          <td>TTGGGCGTCGGTGTCGCGTG</td>
          <td>TTGC</td>
          <td>-0.517138</td>
        </tr>
        <tr>
          <th>3454</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>3454</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TAGATGCTTGGGCATCAATATCACGTGGAACC</td>
          <td>GCGA</td>
          <td>GGAAC</td>
          <td>-14</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>13</td>
          <td>Variant</td>
          <td>CTTGGGCGTCGGTGTCGCGT</td>
          <td>GCGG</td>
          <td>0.002245</td>
        </tr>
      </tbody>
    </table>
    <p>3455 rows × 21 columns</p>
    </div>



Technical replicates show decent LFC correlation.

.. code:: ipython3

    cdata = cdata1 + cdata2


.. code:: ipython3

    cdata





.. parsed-literal::

    Genome Editing Screen comprised of n_guides x n_conditions = 3455 x 12
       guides:    'name', 'Unnamed: 0', 'Target gene/variant', 'Target descriptor', 'Arbitrary number', 'gRNA position category', 'Target base position in gRNA', 'Target base position in reporter', 'BE', 'Group', 'sequence', 'Reporter', 'barcode', '5-nt PAM', 'offset', 'target', 'target_pos', 'Group2', 'masked_sequence', 'masked_barcode', 'bot_top.lfc.median'
       condit:    'index', 'sort', 'replicate'
       condit_m:  
       condit_p:  
       layers:    'edits', 'X_bcmatch'
       uns:       'allele_counts'



You can concatenate different samples with shared guides.

.. code:: ipython3

    cp.concat((cdata1, cdata2))


.. parsed-literal::

    Genome Editing Screen comprised of n_guides x n_conditions = 3455 x 24
       guides:    'name', 'Unnamed: 0', 'Target gene/variant', 'Target descriptor', 'Arbitrary number', 'gRNA position category', 'Target base position in gRNA', 'Target base position in reporter', 'BE', 'Group', 'sequence', 'Reporter', 'barcode', '5-nt PAM', 'offset', 'target', 'target_pos', 'Group2', 'masked_sequence', 'masked_barcode', 'bot_top.lfc.median'
       condit:    'index', 'sort', 'replicate'
       condit_m:  
       condit_p:  
       layers:    'X', 'X_bcmatch', 'edits', 'lognorm_counts', 'lognorm_edits'
       uns:       'allele_counts'



Getting edit rates from allele counts
-------------------------------------

.. code:: ipython3

    cdata.get_edit_rate(normalize_by_editable_base = False,
                       edited_base = "A",
                       editable_base_start = 3,
                       editable_base_end = 8,
                       bcmatch_thres = 10,
                       prior_weight = 1)


.. code:: ipython3

    cdata.uns["edit_counts"] = cdata.get_edit_from_allele()
    


.. code:: ipython3

    cdata.get_edit_mat_from_uns("A", "G", match_target_position = True)
    cdata.get_edit_rate(edited_base = "A", bcmatch_thres = 10)
    plt.hist(cdata.guides.edit_rate, bins=30)
    plt.show()


.. image:: output_34_1.png



Calculating LFC
~~~~~~~~~~~~~~~

.. code:: ipython3

    cdata.log_norm()
    cdata.log_fold_change_aggregate("bot", "top", aggregate_condit = "replicate")

.. code:: ipython3

    cdata.guides




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>name</th>
          <th>Unnamed: 0</th>
          <th>Target gene/variant</th>
          <th>Target descriptor</th>
          <th>Arbitrary number</th>
          <th>gRNA position category</th>
          <th>Target base position in gRNA</th>
          <th>Target base position in reporter</th>
          <th>BE</th>
          <th>Group</th>
          <th>...</th>
          <th>barcode</th>
          <th>5-nt PAM</th>
          <th>offset</th>
          <th>target</th>
          <th>target_pos</th>
          <th>Group2</th>
          <th>masked_sequence</th>
          <th>masked_barcode</th>
          <th>bot_top.lfc.median</th>
          <th>edit_rate</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>CONTROL_1_g1</td>
          <td>0</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g1</td>
          <td>4</td>
          <td>10</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>GTTT</td>
          <td>GGGAG</td>
          <td>-10</td>
          <td>CONTROL_1</td>
          <td>9</td>
          <td>NegCtrl</td>
          <td>CCTGCGCGGTGGGGGGCTTT</td>
          <td>GTTT</td>
          <td>-0.135550</td>
          <td>0.531163</td>
        </tr>
        <tr>
          <th>1</th>
          <td>CONTROL_1_g2</td>
          <td>1</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g2</td>
          <td>5</td>
          <td>11</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>AACA</td>
          <td>TGGGA</td>
          <td>-11</td>
          <td>CONTROL_1</td>
          <td>10</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCTT</td>
          <td>GGCG</td>
          <td>-0.059391</td>
          <td>0.640765</td>
        </tr>
        <tr>
          <th>2</th>
          <td>CONTROL_1_g3</td>
          <td>2</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g3</td>
          <td>5</td>
          <td>12</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>CGCT</td>
          <td>TTGGG</td>
          <td>-12</td>
          <td>CONTROL_1</td>
          <td>11</td>
          <td>NegCtrl</td>
          <td>CCCTGCGCGGTGGGGGGCT</td>
          <td>CGCT</td>
          <td>0.141290</td>
          <td>0.417709</td>
        </tr>
        <tr>
          <th>3</th>
          <td>CONTROL_1_g4</td>
          <td>3</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>TGAG</td>
          <td>TTTGG</td>
          <td>-13</td>
          <td>CONTROL_1</td>
          <td>12</td>
          <td>NegCtrl</td>
          <td>GGCCCTGCGCGGTGGGGGGC</td>
          <td>TGGG</td>
          <td>-0.072358</td>
          <td>0.126400</td>
        </tr>
        <tr>
          <th>4</th>
          <td>CONTROL_1_g5</td>
          <td>4</td>
          <td>CONTROL</td>
          <td>NaN</td>
          <td>1</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>NegCtrl</td>
          <td>...</td>
          <td>GTAT</td>
          <td>CTTTG</td>
          <td>-14</td>
          <td>CONTROL_1</td>
          <td>13</td>
          <td>NegCtrl</td>
          <td>GGGCCCTGCGCGGTGGGGGG</td>
          <td>GTGT</td>
          <td>0.269650</td>
          <td>0.201104</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>3450</th>
          <td>rs9987289_Maj_ABE_347_g1</td>
          <td>3450</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g1</td>
          <td>3</td>
          <td>10</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>CAGT</td>
          <td>CCAGC</td>
          <td>-10</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>9</td>
          <td>Variant</td>
          <td>GCGTCGGTGTCGCGTGGGG</td>
          <td>CGGT</td>
          <td>-0.230264</td>
          <td>0.087379</td>
        </tr>
        <tr>
          <th>3451</th>
          <td>rs9987289_Maj_ABE_347_g2</td>
          <td>3451</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g2</td>
          <td>4</td>
          <td>11</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TCGC</td>
          <td>ACCAG</td>
          <td>-11</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>10</td>
          <td>Variant</td>
          <td>GGCGTCGGTGTCGCGTGGG</td>
          <td>TCGC</td>
          <td>-0.182151</td>
          <td>0.299923</td>
        </tr>
        <tr>
          <th>3452</th>
          <td>rs9987289_Maj_ABE_347_g3</td>
          <td>3452</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g3</td>
          <td>6</td>
          <td>12</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>GCAC</td>
          <td>AACCA</td>
          <td>-12</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>11</td>
          <td>Variant</td>
          <td>TGGGCGTCGGTGTCGCGTGG</td>
          <td>GCGC</td>
          <td>-0.165778</td>
          <td>0.224973</td>
        </tr>
        <tr>
          <th>3453</th>
          <td>rs9987289_Maj_ABE_347_g4</td>
          <td>3453</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g4</td>
          <td>7</td>
          <td>13</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>TTGC</td>
          <td>GAACC</td>
          <td>-13</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>12</td>
          <td>Variant</td>
          <td>TTGGGCGTCGGTGTCGCGTG</td>
          <td>TTGC</td>
          <td>-0.340590</td>
          <td>0.265378</td>
        </tr>
        <tr>
          <th>3454</th>
          <td>rs9987289_Maj_ABE_347_g5</td>
          <td>3454</td>
          <td>rs9987289</td>
          <td>Maj</td>
          <td>347</td>
          <td>g5</td>
          <td>8</td>
          <td>14</td>
          <td>ABE</td>
          <td>Variant</td>
          <td>...</td>
          <td>GCGA</td>
          <td>GGAAC</td>
          <td>-14</td>
          <td>rs9987289_Maj_ABE_347</td>
          <td>13</td>
          <td>Variant</td>
          <td>CTTGGGCGTCGGTGTCGCGT</td>
          <td>GCGG</td>
          <td>0.034365</td>
          <td>0.266573</td>
        </tr>
      </tbody>
    </table>
    <p>3455 rows × 22 columns</p>
    </div>



Allele translation
~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    cdata_tiling = cp.read_h5ad("../../072121_ABE_topbot/beret_counts/LDLRCDS/032422_crispresso/beret_count_072121_ABE_topbot_LDLRCDS.h5ad")


.. code:: ipython3

    cdata_tiling.uns["allele_counts"].allele
    





.. parsed-literal::

    0                                         11224415:14:+:A>G
    1                        11224401:0:+:A>G,11224415:14:+:A>G
    2                        11224410:9:+:A>G,11224415:14:+:A>G
    3         11224401:0:+:A>G,11224402:1:+:A>G,11224410:9:+...
    4                                          11224401:0:+:A>G
                                    ...                        
    438001    11203000:4:+:A>G,11203002:6:+:A>G,11203006:10:...
    438002    11224074:0:+:A>G,11224086:12:+:A>G,11224092:18...
    438003    0:0:+:A>G,3:3:+:A>G,11:11:+:A>G,13:13:+:A>G,17...
    438004                  11217409:23:+:G>-,11217417:31:+:->C
    438005    11226735:30:-:A>G,11226742:23:-:A>G,11226747:1...
    Name: allele, Length: 438006, dtype: object



Writing
~~~~~~~

.. code:: ipython3

    cdata.to_Excel("tmp.xlsx")


.. parsed-literal::

    Writing to: tmp.xlsx
    
    	Sheet 1:	X
    	Sheet 2:	edits
    	Sheet 3:	X_bcmatch
    	Sheet 4:	lognorm_counts
    	Sheet 5:	lognorm_edits
    	Sheet 6:	guides
    	Sheet 7:	condit
    	Sheet 8:	screen.uns.allele_counts
    	Sheet 9:	screen.uns.edit_counts


.. code:: ipython3

    cdata.to_mageck_input("mageck_input.txt", target_column='target')

.. code:: bash

    %%bash
    head mageck_input.txt


.. parsed-literal::

    sgRNA	gene	0	1	2	3	4	5	6	7	8	9	10	11
    CONTROL_1_g1	CONTROL_1	171	451	251	422	573	389	456	420	835	435	794	439
    CONTROL_1_g2	CONTROL_1	145	278	257	206	364	273	389	254	527	498	768	195
    CONTROL_1_g3	CONTROL_1	333	835	488	632	898	899	780	713	1189	626	1146	603
    CONTROL_1_g4	CONTROL_1	246	663	387	448	823	595	705	600	921	595	1143	506
    CONTROL_1_g5	CONTROL_1	243	647	434	529	776	451	700	676	1062	611	928	379
    CONTROL_10_g1	CONTROL_10	138	329	229	213	422	292	432	352	409	243	390	274
    CONTROL_10_g2	CONTROL_10	187	468	402	479	643	369	428	469	796	422	787	404
    CONTROL_10_g3	CONTROL_10	57	126	83	131	281	114	184	115	300	106	299	106
    CONTROL_10_g4	CONTROL_10	66	112	120	136	182	128	169	181	256	144	258	179

