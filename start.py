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
    <tr>
    <td valign="top"><ul>
    <li><a href="{appbase}/isotherm/isotherm.ipynb" target="_blank">Compute one</a>
    <li><a href="{appbase}/isotherm/analyse_results.ipynb" target="_blank">Analyse the results</a>
    </ul></td>
     
     <td valign="top"><ul>
    <li><a href="{appbase}/pores/pores.ipynb" target="_blank">Compute Pores</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/multistage_geo_opt.ipynb" target="_blank">Multistage Geo Opt</a>
    <li><a href="{appbase}/isotherm/isotherm_multi.ipynb" target="_blank">Analyse results</a>
    </ul></td>
    </tr></table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase)
    return ipw.HTML(html)
    
#EOF
