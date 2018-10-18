    def run_geom_zeopp(self):
        """This is the main function that will perform a raspa
        calculation for the current pressure"""
        
        
        # network parameters
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        NetworkParameters = DataFactory('zeopp.parameters')
        params = NetworkParameters({
            'ha': True,
            'res': True,
            'sa': [sigma, sigma, 100000],
            'volpo': [sigma, sigma, 100000],
        }

        # Create the input dictionary
        inputs = {
            'code'       : self.inputs.zeopp_code,
            'structure'  : self.ctx.structure,
            'parameters' : params,
            '_options'   : self.inputs.zeopp_options,
            '_label'     : "run_geom_zeopp",
        }

        # Create the calculation process and launch it
        future  = submit(ZeoppCaclculation.process(), **inputs)
        self.report("pk: {} | Running geometry analysis with zeo++".format(future.pid))

        return ToContext(zeopp=Outputs(future))

    def parse_geom_zeopp(self):
        """
        Extract the pressure and loading average of the last completed raspa calculation
        """
        self.ctx.raspa_parameters['GeneralSettings']['HeliumVoidFraction'] = \
        self.ctx.zeopp["pore_volume_volpo"].dict.POAV_Volume_fraction

    def should_run_block_zeopp(self):
        """If the pore non-accessible volume is 0 - there is no need to run"""
        return self.ctx.zeopp["pore_volume_volpo"].dict.PONAV_Volume_fraction <= 0.001
        
    def run_block_zeopp(self):
        """This is the main function that will perform a raspa
        calculation for the current pressure."""

        # Create the input dictionary
        sigma = self.inputs.probe_molecule.get_dict()['sigma']
        params = NetworkParameters({
            'ha':True,
            'block': [sigma, 200],
        })

        inputs = {
            'code'       : self.inputs.zeopp_code,
            'structure'  : self.ctx.structure,
            'parameters' : params,
            '_options'   : self.inputs._options,
            '_label'     : "run_block_zeopp",

        }

        # Create the calculation process and launch it
        future  = submit(ZeoppCaclulation.process(), **inputs)
        self.report("pk: {} | Running zeo++ block volume calculation".format(future.pid))
        self.ctx.block_pk = future.pid
        return ToContext(zeopp=Outputs(future))
    
    def parse_block_zeopp(self):
        """Extract the pressure and loading average of the last completed raspa calculation"""

        self.ctx.raspa_parameters['Component'][0]['BlockPockets'] = True
        self.ctx.raspa_parameters['Component'][0]['BlockPocketsPk'] = self.ctx.block_pk
