import ipywidgets as ipw

def get_start_widget(appbase, jupbase):
    #http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="text-align:center">Isotherm</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Charges</th>        
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Import your data to AiiDA</th>        
    <tr>
    <td valign="top"><ul>
    <li><a href="{appbase}/isotherm/isotherm.ipynb" target="_blank">Compute one</a>
    <li><a href="{appbase}/isotherm/isotherm_multi.ipynb" target="_blank">Compute multiple</a>
    <li><a href="{appbase}/isotherm/database_analysis.ipynb" target="_blank">Analysis</a>
    </ul></td>
     
     <td valign="top"><ul>
    <li><a href="{appbase}/charges/compute_charges.ipynb" target="_blank">Compute Charges</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/import/import.ipynb" target="_blank">Import database</a>
    <li><a href="{appbase}/import/plot.ipynb" target="_blank">Plot imported data</a>
    </ul></td>
    
    </tr></table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase)
    return ipw.HTML(html)
    
#EOF
