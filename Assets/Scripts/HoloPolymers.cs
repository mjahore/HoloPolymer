using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using UnityEngine;

public class HoloPolymers : MonoBehaviour
{
    // Configuration
    public float boxSize = 30.0f;
    public int polyLen = 20;
    public float initVel = 1.0f;
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
    public float bondLength = 0.7f;
    public float Rc = 2.0f;
    public float Rm = 1.0f;
    public float bondStrength = 100.0f;
    public float bendStrength = 100.0f;
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

    // Private data
    private float scaleFactor;
    private float avgBondLength;
    private float avgBondAngle;
    private float RhConversion = 1.0f / 0.774f;
    private GameObject VisualRg;
    private GameObject VisualRh;
    
    // These are some public routines to allow the user to change the
    // appearance on the fly.
    public void ResetSim()
    {
        for (int i=0; i<polyLen; i++)
        {
            simulationParticles[i].GetComponent<SimParticle>().Reset(scaleFactor);
        }
    }


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


    // Start is called before the first frame update
    void Start()
    {

        // Scale factor?
        scaleFactor = 2.0f * GameObject.Find("BaseParticle").transform.localScale.x / boxSize;

        // System parameters
        float systemVolume = boxSize * boxSize * boxSize;
        sigma = Mathf.Sqrt(6.0f * gamma * temperature / dt);

        // Objects for the scene
        simulationParticles = new GameObject[polyLen];
        bondVectors = new GameObject[polyLen - 1];

        // Set up the polymer chain.
        for (int i = 0; i < polyLen; i++)
        {

            Vector3 initialPosition;
            if (i == 0)
            {
                initialPosition = GameObject.Find("BaseParticle").transform.position;
            }
            else
            {
                //initialPosition = new Vector3(Random.Range(-1, 1), Random.Range(-1, 1), Random.Range(-1, 1));
                initialPosition = new Vector3(1, 0, 0);
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

        // Move the center-of-mass to <0, 0, 0>:
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

    // Note that we use FixedUpdate(), which is framerate independent. This is necessary 
    // to get the correct physics. Default is this is called 50 times per second, so we 
    // can assign a dt of 0.1. This is needed for the random force in DPD.
    void Update()
    {
        // Check for equilibrate flag to speed up
        for (int i = 0; i < polyLen; i++)
        {
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

        if (showBondVectors) UpdateBonds();

        if (showRgSphere)
        {
            UpdateRgSphere();
        } 
        else
        {
            Destroy(VisualRg);
        }

        if (showRhSphere)
        {
            UpdateRhSphere();
        } 
        else
        {
            Destroy(VisualRh);
        }


        FENEBonds();
        if (fixedBondAngles) BondAngles();
        Interparticle_Forces();

        // Correct the velocity according to velocity Verlet algorithm.
        for (int i = 0; i < polyLen; i++)
        {
            simulationParticles[i].GetComponent<SimParticle>().CorrectVelocity();
        }

    }

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


    void UpdateRgSphere()
    {
        if (!VisualRg) VisualRg = Instantiate(rgSphere);
        float Rg = RadiusOfGyration();
        VisualRg.transform.localScale    = new Vector3(Rg, Rg, Rg);
        VisualRg.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
    }


    void UpdateRhSphere()
    {
        if (!VisualRh) VisualRh = Instantiate(rhSphere);
        float Rg = RadiusOfGyration();
        VisualRh.transform.localScale = new Vector3(Rg, Rg, Rg) * RhConversion;
        VisualRh.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
    }


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


    void DrawRgSphere()
    {
        VisualRg = Instantiate(rgSphere);
        VisualRg.transform.localPosition = GameObject.Find("BaseParticle").transform.position;
        float Rg = RadiusOfGyration();
        VisualRg.transform.localScale = new Vector3(Rg, Rg, Rg);

    }


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


    void BondAngles()
    {
        UnityEngine.Vector3 dr0, dr1;
        float dr0_dot_dr1;
        float ct0, ct1;
        float cos_theta;
        UnityEngine.Vector3 f_total;
        float f_coefficient;
        float theta0 = 1.186823891f; // 68 degrees.

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



    void FENEBonds()
    {
        float max_bond2 = Mathf.Pow(3*bondLength, 2);

        avgBondLength = 0.0f;

        for (int i = 0; i < polyLen - 1; i++)
        {

            // Harmonic force = -x(r_i - r_j);
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
        Debug.Log("Bond: " + avgBondLength + "     Angle: " + avgBondAngle + "    Rg: " + RadiusOfGyration());
    }


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
        Debug.Log("Bond: " + avgBondLength + "     Angle: " + avgBondAngle + "    Rg: " + RadiusOfGyration());

    }

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



    /*
     * CenterOfMass() - Returns the center of mass of the polymer chain, and also
     *                  translates the polymer chain to the center of the simulation
     *                  container object. This helps us with moving the hologram.
     */
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