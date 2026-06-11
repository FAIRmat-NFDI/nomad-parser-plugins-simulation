from nomad.parsing.file_parser import Quantity, TextParser


class OutParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity('program_version', r'^LOBSTER *v([\d\.]+) *', repeats=False),
            Quantity(
                'datetime',
                r'starting on host \S* on (\d{4}-\d\d-\d\d\sat\s\d\d:\d\d:\d\d)\s[A-Z]{3,4}',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'x_lobster_code',
                r'detecting used PAW program... (.*)',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'x_lobster_basis',
                r'setting up local basis functions\.\.\.\s*(?:WARNING.*\s*)*\s*((?:[a-zA-Z]{1,2}\s+\(.+\)(?:\s+\d\S+)+\s+)+)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'x_lobster_basis_species',
                            r'([a-zA-Z]+){1,2}\s+\(([^)]+)\)((?:\s+\d\S+)+)\s+',
                            repeats=True,
                        )
                    ]
                ),
            ),
            Quantity(
                'spilling',
                r'((?:spillings|abs. )[\s\S]*?charge\s*spilling:\s*\d+\.\d+%)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'abs_total_spilling',
                            r'abs.\s*total\s*spilling:\s*(\d+\.\d+)%',
                            repeats=False,
                        ),
                        Quantity(
                            'abs_charge_spilling',
                            r'abs.\s*charge\s*spilling:\s*(\d+\.\d+)%',
                            repeats=False,
                        ),
                    ]
                ),
            ),
            Quantity('finished', r'finished in (\d)', repeats=False),
        ]


def get_icoxplist_quantities(version: str) -> list[Quantity]:
    def icoxp_line_split(line) -> list[str | float | int]:
        tmp = line.split()
        # LOBSTER version 3 and above
        if len(tmp) == 8:
            return [
                tmp[0],
                tmp[1],
                tmp[2],
                float(tmp[3]),
                [int(tmp[4]), int(tmp[5]), int(tmp[6])],
                float(tmp[7]),
            ]
        elif len(tmp) == 9 and not tmp[-1].isdigit():
            # Spin polarized data LOBSTER version 5.1 and above
            return [
                tmp[0],
                tmp[1],
                tmp[2],
                float(tmp[3]),
                [int(tmp[4]), int(tmp[5]), int(tmp[6])],
                [float(tmp[7]), float(tmp[8])],
            ]
        elif len(tmp) == 9 and tmp[-1].isdigit():
            # Non-Spin polarized data LOBSTER version 5.1 and above
            return [
                tmp[0],
                tmp[1],
                tmp[2],
                float(tmp[3]),
                [int(tmp[4]), int(tmp[5]), int(tmp[6])],
                float(tmp[7]),
            ]
        # LOBSTER versions below 3
        elif len(tmp) == 6:
            return [tmp[0], tmp[1], tmp[2], float(tmp[3]), float(tmp[4]), int(tmp[5])]

    float_version = float(
        version.split('.', maxsplit=1)[0] + '.' + version.split('.')[1]
    )
    if 5 > float_version > 2:
        return [
            Quantity(
                'icoxpslist_for_spin',
                r' *(CO[O,H,B,I,P]).*spin *\d *([^#]+[-\d\.]+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'line',
                            # LOBSTER version 3 and above
                            r'( *\d+ +[^ ]+ +[^ ]+ +[\.\d]+ +[-\d]+ +[-\d]+ +[-\d]+ +[-\.\d]+ *)',
                            repeats=True,
                            str_operation=icoxp_line_split,
                        )
                    ]
                ),
            )
        ]
    elif float_version >= 5:
        return [
            Quantity(
                'icoxpslist_for_spin',
                r' *(CO[O,H,B,I,P]).*\n .*([^#]+[-\d\.]+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'line',
                            # LOBSTER version 5.1 and above
                            r'( *\d+ +[^ ]+ +[^ ]+ +[\.\d]+ +[-\d]+ +[-\d]+ +[-\d]+ +[-\.\d]+ +[-\.\d]+)',
                            repeats=True,
                            str_operation=icoxp_line_split,
                        )
                    ]
                ),
            ),
            Quantity(
                'icoxpslist_for_nsp',
                r' *(CO[O,H,B,I,P]).*\n .*([^#]+[-\d\.]+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'line',
                            # LOBSTER version 5.1 and above
                            r'( *\d+ +[^ ]+ +[^ ]+ +[\.\d]+ +[-\d]+ +[-\d]+ +[-\d]+ +[-\.\d]+ *)',
                            repeats=True,
                            str_operation=icoxp_line_split,
                        )
                    ]
                ),
            ),
        ]
    else:
        return [
            Quantity(
                'icoxpslist_for_spin',
                r' *(CO[OH]P).*spin *\d *([^#]+[-\d\.]+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'line',
                            # LOBSTER versions below 3
                            r'( *\d+ +[^ ]+ +[^ ]+ +[\.\d]+ +[-\.\d]+ +[\d]+ *)',
                            repeats=True,
                            str_operation=icoxp_line_split,
                        ),
                    ]
                ),
            )
        ]


class ICOXPLISTParser(TextParser):
    version: str

    def init_quantities(self):
        self._quantities = get_icoxplist_quantities(self.version)

    def reset(self):
        super().reset()
        self.init_quantities()
