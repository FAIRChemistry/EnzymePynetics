
    def _initialize_models(
        self, substrate, product, enzyme, inhibitor, only_irrev_MM: bool
    ) -> Dict[str, KineticModel]:
        """Initializer for all kinetic models. If an inhibitor is provided, inhibition models are initialized. If no inhibitor is specified,
        inhibitory models for substrate and product inhibition are initialized additionally to the irreversible Michaelis Menten model.
        """

        y0 = self._get_y0s(
            substrate=substrate, product=product, enzyme=enzyme, inhibitor=inhibitor
        )

        irreversible_Michaelis_Menten_inactivation = KineticModel(
            name="irreversible Michaelis Menten with enzyme inactivation",
            params=[],
            y0=y0,
            kcat_initial=self.initial_kcat,
            Km_initial=self.initial_Km,
            model=irreversible_model,
            enzyme_inactivation=True,
        )

        irreversible_Michaelis_Menten = KineticModel(
            name="irreversible Michaelis Menten",
            params=[],
            y0=y0,
            kcat_initial=self.initial_kcat,
            Km_initial=self.initial_Km,
            model=irreversible_model,
            enzyme_inactivation=False,
        )

        if only_irrev_MM:
            return {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
            }

        if np.all(self.inhibitor == 0):
            # Product inhibition models

            y0 = self._get_y0s(
                substrate=self.substrate,
                enzyme=self.enzyme,
                product=self.product,
                inhibitor=self.product,
            )

            competitive_product_inhibition = KineticModel(
                name="competitive product inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_product_inhibition_inactivation = KineticModel(
                name="competitive product inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_product_inhibition = KineticModel(
                name="uncompetitive product inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_product_inhibition_inactivation = KineticModel(
                name="uncompetitive product inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_product_inhibition = KineticModel(
                name="non-competitive product inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_product_inhibition_inactivation = KineticModel(
                name="non-competitive product inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_product_inhibition_model,
                enzyme_inactivation=True,
            )

            # Substrate inhibition models

            y0 = self._get_y0s(
                substrate=self.substrate,
                enzyme=self.enzyme,
                product=self.product,
                inhibitor=self.substrate,
            )

            substrate_inhibition = KineticModel(
                name="substrate inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            substrate_inhibition_inactivation = KineticModel(
                name="substrate inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=substrate_inhibition_model,
                enzyme_inactivation=True,
            )

            model_dict = {
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
                competitive_product_inhibition.name: competitive_product_inhibition,
                competitive_product_inhibition_inactivation.name: competitive_product_inhibition_inactivation,
                uncompetitive_product_inhibition.name: uncompetitive_product_inhibition,
                uncompetitive_product_inhibition_inactivation.name: uncompetitive_product_inhibition_inactivation,
                uncompetitive_product_inhibition.name: uncompetitive_product_inhibition,
                uncompetitive_product_inhibition_inactivation.name: uncompetitive_product_inhibition_inactivation,
                noncompetitive_product_inhibition.name: noncompetitive_product_inhibition,
                noncompetitive_product_inhibition_inactivation.name: noncompetitive_product_inhibition_inactivation,
                substrate_inhibition.name: substrate_inhibition,
                substrate_inhibition_inactivation.name: substrate_inhibition_inactivation,
            }

            return model_dict

        else:
            y0 = self._get_y0s(
                self.substrate, self.enzyme, self.product, self.inhibitor
            )

            competitive_inhibition = KineticModel(
                name="competitive inhibition",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            competitive_inhibition_inactivation = KineticModel(
                name="competitive inhibition with enzyme inactivation",
                params=["K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_model,
                enzyme_inactivation=True,
            )

            uncompetitive_inhibition = KineticModel(
                name="uncompetitive inhibition",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            uncompetitive_inhibition_inactivation = KineticModel(
                name="uncompetitive inhibition with enzyme inactivation",
                params=["K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=uncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            noncompetitive_inhibition = KineticModel(
                name="non-competitive inhibition",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=False,
            )
            noncompetitive_inhibition_inactivation = KineticModel(
                name="non-competitive inhibition with enzyme inactivation",
                params=["K_iu", "K_ic"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=noncompetitive_inhibition_model,
                enzyme_inactivation=True,
            )

            partially_competitive_inhibition = KineticModel(
                name="partially competitive inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=False,
            )
            partially_competitive_inhibition_inactivation = KineticModel(
                name="partially competitive inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=partially_competitive_inhibition_model,
                enzyme_inactivation=True,
            )
            competitive_inhibition_with_substrate_inhibition_model

            comp_inhib_sub_inhib = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_with_substrate_inhibition_model,
                enzyme_inactivation=False,
            )
            comp_inhib_sub_inhib_enz_inact = KineticModel(
                name="competitive inhibition with uncompetitive substrate inhibition with enzyme inactivation",
                params=["K_ic", "K_iu"],
                y0=y0,
                kcat_initial=self.initial_kcat,
                Km_initial=self.initial_Km,
                model=competitive_inhibition_with_substrate_inhibition_model,
                enzyme_inactivation=True,
            )

            return {
                comp_inhib_sub_inhib_enz_inact.name: comp_inhib_sub_inhib_enz_inact,
                comp_inhib_sub_inhib.name: comp_inhib_sub_inhib,
                irreversible_Michaelis_Menten.name: irreversible_Michaelis_Menten,
                irreversible_Michaelis_Menten_inactivation.name: irreversible_Michaelis_Menten_inactivation,
                competitive_inhibition.name: competitive_inhibition,
                competitive_inhibition_inactivation.name: competitive_inhibition_inactivation,
                uncompetitive_inhibition.name: uncompetitive_inhibition,
                uncompetitive_inhibition_inactivation.name: uncompetitive_inhibition_inactivation,
                # noncompetitive_inhibition.name: noncompetitive_inhibition,
                # noncompetitive_inhibition_inactivation.name: noncompetitive_inhibition_inactivation,
                partially_competitive_inhibition.name: partially_competitive_inhibition,
                partially_competitive_inhibition_inactivation.name: partially_competitive_inhibition_inactivation,
            }
