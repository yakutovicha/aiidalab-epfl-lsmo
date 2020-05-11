# Setting up AiiDA computers & codes

This folder contains a number of [YAML](https://yaml.org/) configuration files in order to set up AiiDA computers and codes,
which can be used automatically with the top-level `setup.py`.

Alternatively, you can use them to set up computers one by one:
```
verdi computer setup --config fidis/fidis.yml
verdi code setup --config cp2k-5.1@fidis.yml
```

See also the AiiDA docs on [setting up a new computer](https://aiida-core.readthedocs.io/en/latest/get_started/computers.html)
and  [setting up a new code](https://aiida-core.readthedocs.io/en/latest/get_started/codes.html).

## For maintainers: create export file with codes & computers

Create export file for others members to import:
```
python create_export.py
```

Note: This uses a temporary AiiDA profile and won't pollute your database.
