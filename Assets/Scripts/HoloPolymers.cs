/*
 * HoloPolymers.cs - An implementation of an implicit solvent simulation of a single polymer chain
 *                   in Unity. The model used here is due to M. Laradji et al. and can be found in 
 *                   the following publication:
 *
 *                   J.D. Revalee, M. Laradji, and P. B. S. Kumar, J. Chem. Phys. 2008, 128, 01B614.
 *
 *                   The thermostat is a Langevin thermostat.
 *
 * Mike Hore (hore@case.edu), Case Western Reserve University
 * July 2020
 *
 *
 * Copyright 2020 Michael J. A. Hore
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using UnityEngine;

public class HoloPolymers : MonoBehaviour
{
    // Configuration
    public float boxSize = 50.0f;
    public int initialPolyLen = 50;
    public float initVel = 1.5f;
    public bool showBondVectors = true;
    public bool showRgSphere = true;
    public bool showRhSphere = true;
    public bool fixedBondAngles = false;
    public bool equilibrate = true;

    // Prefabs
    public GameObject monomerObject;
    public GameObject rgSphere;
    public GameObject rhSphere;
    public GameObject bondVector;

    // Interaction parameters
    public float bondLength = 1.5f;
    public float Rc = 2.0f;
    public float Rm = 1.0f;
    public float bondStrength = 100.0f;
    public float bendStrength = 50.0f;
    public float U_MAX = 200.0f;
    public float U_MIN = -5.0f;

    // Langevin thermostat
    public float temperature = 3.0f;
    public float gamma = 1.0f;
    public float sigma;
    public float dt = 0.0005f;

    // Private globals
    private GameObject[] simulationParticles;
    private GameObject[] bondVectors;
    private float scaleFactor;
    private float avgBondLength;
    private float avgBondAngle;
    private int   polyLen;
    private float RhConversion = 1.2919896f; // = 1.0/0.774
    private int MIN_MONOMERS = 5;
    private int MAX_MONOMERS = 200;
    private GameObject VisualRg;
    private GameObject VisualRh;
    
    ////
    // These are some public routines to allow the user to change the
    // appearance on the fly.
    ////

    //
    // ResetSim() - Returns all monomers to their initial positions/velocities.
    //
    public void ResetSim()
    {
        for (int i=0; i<polyLen; i++)
        {
            simulationParticles[i].GetComponent<SimParticle>().Reset(scaleFactor);
        }
    }


    //
    // ToggleAngles() - Toggles between the freely jointed chain (FJC) model, where any bond angle is 
    //                  permitted, and the freely rotating chain (FRC) model, where a fixed bond angle
    //                  is enforced through a bending potential. By default, this angle is 68 deg. in
    //                  this code (hard coded below). 
    public void ToggleAngles()
    {
        if (fixedBondAngles)
        {
            fixedBondAngles = false;
        }
        else
        {
            fixedBondAngles = true;
        }
    }


    //
    // ToggleBondVectors() - Toggles between displaying and not displaying the bonds between monomers.
    //
    public void ToggleBondVectors()
    {
        if (showBondVectors)
        {
            showBondVectors = false;
            for (int i=0; i<polyLen - 1; i++)
            {
                bondVectors[i].GetComponent<Renderer>().enabled = false;
            }
        } 
        else
        {
            showBondVectors = true;
            for (int i = 0; i < polyLen - 1; i++)
            {
                bondVectors[i].GetComponent<Renderer>().enabled = true;
            }
        } 

    }


    //
    // ToggleRh() - Toggles between displaying and not displaying a sphere representing the
    //              hydrodynamic radius of the polymer chain.
    //
    public void ToggleRh()
    {
        if (showRhSphere)
        {
            showRhSphere = false;
        }
        else
        {
            showRhSphere = true;
        }
    }


    //
    // ToggleRg() - Toggles between displaying and not displaying a sphere representing the
    //              radius of gyration of the polymer chain.
    //
    public void ToggleRg()
    {
        if (showRgSphere)
        {
            showRgSphere = false;
        }
        else
        {
            showRgSphere = true;
        }
    }


    //
    // ToggleSpeed() - If the timestep dt is too large, the simulation will blow up. This is
    //                 especially true when the FRC model is used. ToggleSpeed() switches between
    //                 a timestep of dt, or 10*dt. In the FRJ chain model, using the default parameters,
    //                 a timestep of up to dt = 0.01 is usually stable.
    //
    public void ToggleSpeed()
    {
        if (equilibrate)
        {
            equilibrate = false;
        } 
        else
        {
            equilibrate = true;
        }
    }


    //
    // ScalePolyLen() - Adds or removes monomers from the system when the bounding box object scales
    //                  the BaseParticle object. This will change polyLen. We leave scaleFactor as is.
    //
    public void ScalePolyLen()
    {
        GameObject baseParticle = GameObject.Find("BaseParticle");

        // Calculate the number of monomers based on the size of the bounding box. If localScale = 1,
        // then the box is at the initial size, so we should recover initialPolyLen;
        int newPolyLen = Mathf.RoundToInt(initialPolyLen * baseParticle.transform.localScale.x/0.5f);

        // Safety checks!
        if (newPolyLen > MAX_MONOMERS)
        {
            newPolyLen = MAX_MONOMERS;
        } 

        if (newPolyLen < MIN_MONOMERS)
        {
            newPolyLen = MIN_MONOMERS;
        }


        // Either remove or add monomers to the system.
        if (newPolyLen < polyLen)
        {
            // Remove:
            for (int i = polyLen-1; i > newPolyLen-1; i--)
            {
                // Remove the monomer from the simulation.
                Destroy(simulationParticles[i]);

                // Remove the bond vector from the simulation.
                Destroy(bondVectors[i - 1]);
            }
        }
        else
        {
            for (int i = polyLen; i < newPolyLen; i++)
            {

                // Add:
                Vector3 initialPosition = new Vector3(Random.Range(-1, 1), Random.Range(-1, 1), Random.Range(-1, 1));
                initialPosition = simulationParticles[i - 1].GetComponent<SimParticle>().realPosition + bondLength * initialPosition.normalized;


                // Add the bead, set initial position.
                simulationParticles[i] = Instantiate(monomerObject);
                simulationParticles[i].GetComponent<SimParticle>().MoveBead(initialPosition, scaleFactor);

                // Assign an id.
                simulationParticles[i].GetComponent<SimParticle>().id = i;

                // Species/particle type
                simulationParticles[i].GetComponent<SimParticle>().species = 0;

                // Periodic boundary condition
                simulationParticles[i].GetComponent<SimParticle>().SetBoxSize(boxSize);

                // Set dt
                simulationParticles[i].GetComponent<SimParticle>().SetDt(dt);

                // Set the initial velocity to be that of the previous monomer to keep things near
                // equilibrium.
                simulationParticles[i].GetComponent<SimParticle>().SetVelocity(simulationParticles[i - 1].GetComponent<SimParticle>().beadVelocity);

                // Make small!
                simulationParticles[i].transform.localScale = new Vector3(scaleFactor, scaleFactor, scaleFactor);

                // Add the bond vectors
                UnityEngine.Vector3 offset = simulationParticles[i].transform.localPosition - simulationParticles[i-1].transform.localPosition;
                UnityEngine.Vector3 scale = new Vector3(scaleFactor * 0.1f, offset.magnitude / 2.0f, scaleFactor * 0.1f);
                UnityEngine.Vector3 position = simulationParticles[i-1].transform.localPosition + (offset / 2.0f);
                bondVectors[i-1] = Instantiate(bondVector, position, Quaternion.identity);
                bondVectors[i-1].transform.up = offset;
                bondVectors[i-1].transform.localScale = scale;
            }
        }

        polyLen = newPolyLen;
    }



    
    //
    // Start() - This method is called first, and sets up the simulation. Notice thaat the lengthscales here are
    //           rescaled by scaleFactor when the particles need to be displayed to the user. Thus,they're really
    //           tied to the size of the BoundingBox around the object "BaseParticle", which is the only object
    //           placed in the scene (other than camera, light, and context menu).
    //
    void Start()
    {
        // Initial chain length
        polyLen = initialPolyLen;

        // Scale factor?
        scaleFactor = 2.0f * GameObject.Find("BaseParticle").transform.localScale.x / boxSize;

        // System parameters
        float systemVolume = boxSize * boxSize * boxSize;
        sigma = Mathf.Sqrt(6.0f * gamma * temperature / dt);

        // Objects for the scene. Notice that we allow for a polymer chain having
        // up to MAX_MONOMERS monomers. The number of monomers can be changed on the fly
        // by the user by rescaling the BoundingBox around "BaseParticle".
        simulationParticles = new GameObject[MAX_MONOMERS];
        bondVectors = new GameObject[MAX_MONOMERS-1];

        // Set up the polymer chain, by either anchoring the first monomer to the position of
        // the BaseParticle, or by randomly placing a monomer near the previous monomer in the
        // chain.
        for (int i = 0; i < polyLen; i++)
        {
            Vector3 initialPosition;
            if (i == 0)
            {
                initialPosition = GameObject.Find("BaseParticle").transform.position;
            }
            else
            {
                initialPosition = new Vector3(Random.Range(-1, 1), Random.Range(-1, 1), Random.Range(-1, 1));
                initialPosition = simulationParticles[i - 1].GetComponent<SimParticle>().realPosition + bondLength * initialPosition.normalized;
            }

            // Add the bead, set initial position.
            simulationParticles[i] = Instantiate(monomerObject);
            simulationParticles[i].GetComponent<SimParticle>().MoveBead(initialPosition, scaleFactor);

            // Assign an id.
            simulationParticles[i].GetComponent<SimParticle>().id = i;

            // Species/particle type
            simulationParticles[i].GetComponent<SimParticle>().species = 0;

            // Periodic boundary condition
            simulationParticles[i].GetComponent<SimParticle>().SetBoxSize(boxSize);

            // Set dt
            simulationParticles[i].GetComponent<SimParticle>().SetDt(dt);

            // Get a random velocity for the particle.
            Vector3 initVelocityVector = new Vector3(Random.Range(-1, 1), Random.Range(-1, 1), Random.Range(-1, 1));
            simulationParticles[i].GetComponent<SimParticle>().SetVelocity(initVel * initVelocityVector.normalized);

            // Make small!
            simulationParticles[i].transform.localScale = new Vector3(scaleFactor, scaleFactor, scaleFactor);
        }

        // Move the center-of-mass to <0, 0, 0> to center the chain in the BoundingBox.
        CenterOfMass();

        // Calculate initial forces.
        FENEBonds();
        if (fixedBondAngles) BondAngles();
        Interparticle_Forces();

        // Hide base particle.
        GameObject.Find("BaseParticle").GetComponent<Renderer>().enabled = false;

        // Display bond vectors, Rg, Rh if enabled.
        if (showBondVectors) DrawBonds();
        if (showRgSphere) DrawRgSphere();
        if (showRhSphere) DrawRhSphere();

    }


    //
    // Update() - This is the main loop of the simulation that updates the positions/velocities of the 
    //            monomers using the values of the forces that are calculated. Integration is done using
    //            the velocity-Verlet (VV) algorithm.
    //
    void Update()
    {

	// This is the first step of the VV algorithm where we move the particles, predict a value
        // for the velocity, and apply the Langevin thermostat.
        for (int i = 0; i < polyLen; i++)
        {
            // Check for equilibrate flag to speed up
            if (equilibrate)
            {
                simulationParticles[i].GetComponent<SimParticle>().SetDt(10*dt);
            } else
            {
                simulationParticles[i].GetComponent<SimParticle>().SetDt(dt);
            }
            simulationParticles[i].GetComponent<SimParticle>().MoveBead(scaleFactor);
            simulationParticles[i].GetComponent<SimParticle>().PredictVelocity();
            simulationParticles[i].GetComponent<SimParticle>().Langevin(sigma, gamma);
        }

        // Keep the polymer centered in the scene.
        CenterOfMass();

        // This updates the bond vectors between monomers.
        if (showBondVectors) UpdateBonds();

        // This just shows (or doesn't) a sphere representing Rg.
        if (showRgSphere)
        {
            UpdateRgSphere();
        } 
        else
        {
            Destroy(VisualRg);
        }

        // This is the same, but for Rh.
        if (showRhSphere)
        {
            UpdateRhSphere();
        } 
        else
        {
            Destroy(VisualRh);
        }

        // Recalculate the new values for the force on each bead using the new
        // position of the particles.
        FENEBonds();
        if (fixedBondAngles) BondAngles();
        Interparticle_Forces();

        // Correct the velocity according to velocity Verlet algorithm.
        for (int i = 0; i < polyLen; i++)
        {
            simulationParticles[i].GetComponent<SimParticle>().CorrectVelocity();
        }

    }


    // 
    // UpdateBonds() - This updates the position of the bonds between monomers. A bond vector is just represented by
    //                 a cylinder prefab. As with positions, everything is scaled by scaleFactor.
    //
    void UpdateBonds()
    {
        for (int i = 0; i < polyLen - 1; i++)
        {
            UnityEngine.Vector3 offset = simulationParticles[i + 1].transform.localPosition - simulationParticles[i].transform.localPosition;
            UnityEngine.Vector3 scale = new Vector3(scaleFactor * 0.1f, offset.magnitude / 2.0f, scaleFactor * 0.1f);
            UnityEngine.Vector3 position = simulationParticles[i].transform.localPosition + (offset / 2.0f);
            bondVectors[i].transform.localPosition = position;
            bondVectors[i].transform.up = offset;
            bondVectors[i].transform.localScale = scale;
        }
    }


    //
    // UpdateRgSphere() - As the polymer chain's conformation evolves, its radius of gyration changes. This
    //                    updates the size of the Rg sphere as Rg changes.
    //
    void UpdateRgSphere()
    {
        if (!VisualRg) VisualRg = Instantiate(rgSphere);
        float Rg = RadiusOfGyration();
        VisualRg.transform.localScale    = new Vector3(Rg, Rg, Rg);
        VisualRg.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
    }


    //
    // UpdateRgSphere() - As the polymer chain's conformation evolves, its hydrodynamic radius changes. This
    //                    updates the size of the Rh sphere as Rh changes.
    //
    void UpdateRhSphere()
    {
        if (!VisualRh) VisualRh = Instantiate(rhSphere);

	// Note how we calculate Rh. We are assuming that Rg/Rh = 0.774.
        float Rg = RadiusOfGyration();
        VisualRh.transform.localScale = new Vector3(Rg, Rg, Rg) * RhConversion;

        VisualRh.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
    }


    //
    // DrawBonds() - This method actually creates the bond vectors between monomers, and sets them to their
    //               initial positions & orientations.
    //
    void DrawBonds()
    {
        for (int i = 0; i < polyLen - 1; i++)
        {
            UnityEngine.Vector3 offset = simulationParticles[i + 1].transform.localPosition - simulationParticles[i].transform.localPosition;
            UnityEngine.Vector3 scale = new Vector3(scaleFactor * 0.1f, offset.magnitude / 2.0f, scaleFactor * 0.1f);
            UnityEngine.Vector3 position = simulationParticles[i].transform.localPosition + (offset / 2.0f);
            bondVectors[i] = Instantiate(bondVector, position, Quaternion.identity);
            bondVectors[i].transform.up = offset;
            bondVectors[i].transform.localScale = scale;
        }
    }


    // 
    // DrawRgSphere() - This is a lot like UpdateRgSphere, except that it creates the RgSphere object.
    //
    void DrawRgSphere()
    {
        VisualRg = Instantiate(rgSphere);
        VisualRg.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
        float Rg = RadiusOfGyration();
        VisualRg.transform.localScale = new Vector3(Rg, Rg, Rg);
    }

   
    // 
    // DrawRhSphere() - This is a lot like UpdateRhSphere, except that it creates the RhSphere object.
    //
    void DrawRhSphere()
    {
        if (!VisualRh) VisualRh = Instantiate(rhSphere);
        VisualRh = Instantiate(rhSphere);
        VisualRh.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
        float Rg = RadiusOfGyration();
        VisualRh.transform.localScale = new Vector3(Rg, Rg, Rg) * RhConversion;
    }


    /***********************************************************************
     * Below here is where the physics happens!
     ***********************************************************************/


    //
    // BondAngles() - This applies a potential U = -k * (cos(theta) - cos(theta0))^2 to each
    //                monomer. It's just a derivative and a lot of dot products. If you don't
    //                understand the math, e-mail me.
    void BondAngles()
    {
        UnityEngine.Vector3 dr0, dr1;
        float dr0_dot_dr1;
        float ct0, ct1;
        float cos_theta;
        UnityEngine.Vector3 f_total;
        float f_coefficient;
        float theta0 = 1.186823891f; // 68 degrees in radians.

        avgBondAngle = 0.0f;
        for (int i = 0; i < polyLen - 2; i++)
        {
            dr0 = simulationParticles[i].transform.localPosition - simulationParticles[i + 1].transform.localPosition;
            dr1 = simulationParticles[i + 1].transform.localPosition - simulationParticles[i + 2].transform.localPosition;

            // This gets complicated!
            dr0_dot_dr1 = UnityEngine.Vector3.Dot(dr0, dr1);

            ct0 = dr0_dot_dr1 / Mathf.Pow(dr0.magnitude, 2);
            ct1 = dr0_dot_dr1 / Mathf.Pow(dr1.magnitude, 2);

            cos_theta = dr0_dot_dr1 / dr0.magnitude / dr1.magnitude;

            f_coefficient = bendStrength * (Mathf.Cos(theta0) - cos_theta) / (dr0.magnitude * dr1.magnitude);
            avgBondAngle += cos_theta;

            // First particle of the three-body set:
            f_total = f_coefficient * (dr1 - ct0 * dr0);
            simulationParticles[i].GetComponent<SimParticle>().AddForce(f_total);

            // Second particle of the three-body set:
            f_total = f_coefficient * (ct0 * dr0 - ct1 * dr1 + dr0 - dr1);
            simulationParticles[i + 1].GetComponent<SimParticle>().AddForce(f_total);

            // Third particle of the three-body set:
            f_total = f_coefficient * (ct1 * dr1 - dr0);
            simulationParticles[i + 2].GetComponent<SimParticle>().AddForce(f_total);
        }

        // Convert to degrees:
        avgBondAngle /= (polyLen - 2);
        avgBondAngle = Mathf.Acos(avgBondAngle) * (180.0f/3.14159f);
    }

 
    //
    // FENEBonds() - This applies the FENE potential to each monomer. The FENE potential is a steeper
    //               bonding potential than the harmonic potential, so it enforces the bond lengths better
    //               in my experience.
    //
    void FENEBonds()
    {
        float max_bond2 = Mathf.Pow(3*bondLength, 2);

        avgBondLength = 0.0f;

        for (int i = 0; i < polyLen - 1; i++)
        {

            // FENE force = -grad(FENE potential)
            Vector3 bondVector = simulationParticles[i].GetComponent<SimParticle>().realPosition - simulationParticles[i + 1].GetComponent<SimParticle>().realPosition;
            float harmonicForce = -bondStrength / (1.0f - Mathf.Pow(bondVector.magnitude - bondLength, 2) / max_bond2) * (bondVector.magnitude - bondLength);

            // Add the force to particle #1:
            simulationParticles[i].GetComponent<SimParticle>().AddForce(harmonicForce * bondVector.normalized);

            // Newton's 3rd Law! 
            simulationParticles[i + 1].GetComponent<SimParticle>().AddForce(-harmonicForce * bondVector.normalized);

            // Average bond length:
            avgBondLength += bondVector.magnitude;
        }

        avgBondLength /= (polyLen - 1);
        avgBondLength *= scaleFactor;
        Debug.Log("Monomers:" + polyLen + "     Bond: " + avgBondLength + "     Angle: " + avgBondAngle + "    Rg: " + RadiusOfGyration());
    }


    //
    // HarmonicBonds() - This is another bonding approach, but uses a harmonic potential instead. It's here,
    //                   but I don't recommend using it unless you have a good reason.
    //
    void HarmonicBonds()
    {
        avgBondLength = 0.0f;

        for (int i = 0; i < polyLen - 1; i++)
        {

            // Harmonic force = -x(r_i - r_j);
            Vector3 bondVector = simulationParticles[i].GetComponent<SimParticle>().realPosition - simulationParticles[i + 1].GetComponent<SimParticle>().realPosition;
            float harmonicForce = -bondStrength * (bondVector.magnitude - bondLength);

            // Add the force to particle #1:
            simulationParticles[i].GetComponent<SimParticle>().AddForce(harmonicForce * bondVector.normalized);

            // Newton's 3rd Law! 
            simulationParticles[i + 1].GetComponent<SimParticle>().AddForce(-harmonicForce * bondVector.normalized);

            // Average bond length:
            avgBondLength += bondVector.magnitude;
        }

        avgBondLength /= (polyLen - 1);
        avgBondLength *= scaleFactor;
        Debug.Log("Monomers:" + polyLen + "     Bond: " + avgBondLength + "     Angle: " + avgBondAngle + "    Rg: " + RadiusOfGyration());
    }


    //
    // InterParticle_Forces() - This calculates the forces between all beads according to the model from Laradji.
    //                          Note that I am just brute forcing this. If you had more particles, it would be
    //                          ridiculously slow to do this double loop. You should use a neighbor list instead.
    //
    void Interparticle_Forces()
    {
        int i, j;
        Vector3 f_total, dr, xform;
        float RcRm = Rc - Rm;

        f_total = new Vector3(0, 0, 0);

        for (i = 0; i < polyLen; i++)
        {
            for (j=0; j<polyLen; j++)
            {
                if (j != i)
                {
                    dr = simulationParticles[i].GetComponent<SimParticle>().realPosition - simulationParticles[j].GetComponent<SimParticle>().realPosition;

                    // Check periodic boundary conditions:
                    xform = new Vector3(0, 0, 0);

                    if (dr.magnitude < Rc)
                    {
                        if (dr.magnitude < Rm)
                        {
                            f_total = 2.0f * (U_MAX - U_MIN) * (Rm - dr.magnitude) / Mathf.Pow(Rm, 2) * dr.normalized;
                        }
                        else
                        {
                            f_total = (-6.0f * U_MIN * Mathf.Pow(Rc - dr.magnitude, 2) / Mathf.Pow(RcRm, 3) + 6.0f * U_MIN * (Rc - dr.magnitude) / Mathf.Pow(RcRm, 2)) * dr.normalized;
                        }
                    }
                    else
                    {
                        f_total.Set(0.0f, 0.0f, 0.0f);
                    }

                    simulationParticles[i].GetComponent<SimParticle>().AddForce(f_total);
                }
            }

        }
    }


    //
    // CenterOfMass() - Returns the center of mass of the polymer chain, and also
    //                  translates the polymer chain to the center of the simulation
    //                  container object. This helps us with moving the hologram.
    //
    Vector3 CenterOfMass()
    {
        Vector3 com = new Vector3(0, 0, 0);

        // First calculate center of mass of the polymer chain.
        for (int i = 0; i<polyLen; i++)
        {
            com += simulationParticles[i].transform.localPosition;
        }
        com /= polyLen;

        // We want to translate the polymer back to the center of the empty
        // simulation container object.
        com = com - GameObject.Find("BaseParticle").transform.position;

        // Translate the polymer back to the center.
        for (int i=0; i<polyLen; i++)
        {
            simulationParticles[i].transform.localPosition -= com;
        }
        return com;
    }


    //
    // RadiusOfGyration() - Calculates Rg for the chain.
    //
    float RadiusOfGyration()
    {
        int i, j;
        float Rg;
        Vector3 dr = new Vector3(0, 0, 0);

        Rg = 0.0f;

        for (i = 0; i < polyLen; i++)
        {
            for (j = 0; j < polyLen; j++)
            {
                dr = simulationParticles[i].transform.localPosition - simulationParticles[j].transform.localPosition;
                Rg += Mathf.Pow(dr.magnitude, 2);
            }
        }

        Rg = Mathf.Sqrt(Rg / polyLen / polyLen);

        return Rg;
    }

}
