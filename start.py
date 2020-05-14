import ipywidgets as ipw

def get_start_widget(appbase, jupbase):
    #http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="text-align:center">Isotherm</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Pore analysis</th>        
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Geometry Optimization</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Electronic Properties</th>
    <tr>
    <td valign="top"><ul>
    <li><a href="{appbase}/isotherm/isotherm.ipynb" target="_blank">Compute one</a>
    <li><a href="{appbase}/isotherm/analyse_results.ipynb" target="_blank">Analyse the results</a>
    </ul></td>
     
    <td valign="top"><ul>
    <li><a href="{appbase}/pores/pore_analysis.ipynb" target="_blank">Pore Analysis</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/multistage_geo_opt.ipynb" target="_blank">Geometry Optimization</a>
    <li><a href="{appbase}/multistage_geo_opt_ddec.ipynb" target="_blank">Geometry Optimization and Charges</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/electronic/band_structure.ipynb" target="_blank">Band Structure</a>
    </ul></td>
    
    </tr></table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase)
    return ipw.HTML(html)
    
#EOF
