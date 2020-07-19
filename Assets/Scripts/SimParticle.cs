/*
 * SimParticle.cs - This class adds simulation functionality to a Monomer prefab.
 *
 * Mike Hore (hore@case.edu), Case Western Reserve University
 * July 2020
 *
 */
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Diagnostics;
using System.Security.Cryptography;
using UnityEngine;

public class SimParticle : MonoBehaviour
{
    // Store velocity and forces for current timestep.
    public  Vector3 beadVelocity;
    private Vector3 beadNetForce;

    // Store velocity and forces for previous timestep.
    private Vector3 prevVelocity;
    private Vector3 prevNetForce;

    // Initial coordinates so we can reset the simulation.
    private Vector3 initialPosition;
    private Vector3 initialVelocity;

    // Integration timestep.
    private float dt;

    // Box size for periodic boundary condition.
    private float boxSize;

    // Unique identifier for each particle.
    public int id;

    // Species of this particle. 0 = polymer, 1 = solvent.
    public int species;

    // Coordinates not subject to periodic boundary conditions.
    public Vector3 realPosition;

    // Start is called before the first frame update
    void Start()
    {
        // Initialize. We'll set these in HoloDPD (because
        // it's easier to create a polymer there).
        this.beadVelocity.Set(0.0f, 0.0f, 0.0f);
        this.beadNetForce.Set(0.0f, 0.0f, 0.0f);
        this.prevVelocity.Set(0.0f, 0.0f, 0.0f);
        this.prevNetForce.Set(0.0f, 0.0f, 0.0f);
    }

    // Adds a force vector to this particle's net force.   
    public void AddForce(Vector3 myForce)
    {
        this.beadNetForce += myForce;
    }

    // Correct the velocity based on current force values.
    public void CorrectVelocity()
    {
        this.beadVelocity = this.prevVelocity + 0.5f * dt * (this.beadNetForce + this.prevNetForce);
    }

    // This applies the Langevin thermostat, and initializes the net force vector before the
    // forces are calculated.
    public void Langevin(float sigma, float gamma)
    {
        // Store current velocity/net force.
        this.prevVelocity = this.beadVelocity;
        this.prevNetForce = this.beadNetForce;

        // Reset current net force.
        Vector3 randF = new Vector3(Random.Range(-0.5f, 0.5f), Random.Range(-0.5f, 0.5f), Random.Range(-0.5f, 0.5f));
        this.beadNetForce = (2.0f * sigma * randF - gamma * this.beadVelocity);
    }


    // Once everything's calculated, update position of the bead
    // based on the current velocity and forces.
    public void MoveBead(float scaleFactor)
    {
        // Update the bead position.
        Vector3 dr = (dt * beadVelocity) + (0.5f * dt * dt * beadNetForce);
        this.realPosition = this.realPosition + dr;

        // For display -- we gotta rescale!
        this.transform.localPosition = this.realPosition * scaleFactor;

    }

    // Move the bead to the specified location. Assumption is that this is only called
    // when the monomer object is created. So initial position is defined based on this
    // function.
    public void MoveBead(Vector3 myPosition, float scaleFactor)
    {
        this.GetComponent<Renderer>().transform.localPosition = myPosition * scaleFactor;
        this.realPosition = myPosition;
        this.initialPosition = myPosition;
    }

    // Reset puts everything back to where it was when the simulation was started.
    public void Reset(float scaleFactor)
    {
        // Re-initialize
        this.beadNetForce.Set(0.0f, 0.0f, 0.0f);
        this.prevNetForce.Set(0.0f, 0.0f, 0.0f);

        // Reset the velocity.
        this.beadVelocity = initialVelocity;

        // Move us back to the start!
        this.realPosition = initialPosition;

        // For display.
        this.transform.localPosition = this.realPosition * scaleFactor;
    }


    // Predict the velocity based on current force values.
    public void PredictVelocity()
    {
        Vector3 dv = this.beadVelocity + 0.5f * dt * this.beadNetForce;
        this.beadVelocity = dv;
    }

    // Set box size for periodic boundary conditions.
    public void SetBoxSize (float myBoxSize)
    {
        this.boxSize = myBoxSize;
    }

    // Set integration timestep.
    public void SetDt(float mydt)
    {
        this.dt = mydt;
    }

    // Set the initial velocity.
    public void SetVelocity(Vector3 myVelocity)
    {
        this.beadVelocity = myVelocity;
        this.initialVelocity = myVelocity;
    }

}
