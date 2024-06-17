using UnityEngine;
using System.Collections.Generic;

public class RRT : MonoBehaviour
{
    //public Text statusText; // UI text to display status (optional)
    public GameObject startWaypoint; // Start waypoint assigned in the Inspector
    public GameObject endWaypoint; // End waypoint assigned in the Inspector
    public List<GameObject> obstacles; // List of all obstacle objects assigned in the Inspector
    private Dictionary<GameObject, Dictionary<GameObject, float>> graph = new Dictionary<GameObject, Dictionary<GameObject, float>>();

    void Start()
    {
        Debug.Log("Start: Checking assignments before initialization");
        Debug.Log("startWaypoint: " + (startWaypoint != null ? startWaypoint.name : "null"));
        Debug.Log("endWaypoint: " + (endWaypoint != null ? endWaypoint.name : "null"));

        AddMeshCollidersToObstacles();
        AddMeshCollidersToPoints();

        if (startWaypoint == null || endWaypoint == null)
        {
            Debug.LogError("Start or End waypoint is not assigned. Please assign waypoints in the Inspector.");
            return;
        }

        InitializeGraph();
        Debug.Log("Graph initialized.");

        // Automatically find the path at the start
        FindPath(startWaypoint, endWaypoint);
    }

    void AddMeshCollidersToObstacles()
    {
        foreach (var obstacle in obstacles)
        {
            MeshCollider meshCollider = obstacle.GetComponent<MeshCollider>();
            if (meshCollider == null)
            {
                meshCollider = obstacle.AddComponent<MeshCollider>();
                Debug.Log("Mesh Collider added to " + obstacle.name);
            }
        }
    }

    void AddMeshCollidersToPoints()
    {
        List<GameObject> points = new List<GameObject> { startWaypoint, endWaypoint };

        foreach (var point in points)
        {
            if (point != null)
            {
                MeshCollider meshCollider = point.GetComponent<MeshCollider>();
                if (meshCollider == null)
                {
                    meshCollider = point.AddComponent<MeshCollider>();
                    Debug.Log("Mesh Collider added to " + point.name);
                }
            }
            else
            {
                Debug.LogWarning("Point is null.");
            }
        }
    }

    void InitializeGraph()
    {
        Debug.Log("Initializing graph...");
        graph[startWaypoint] = new Dictionary<GameObject, float>();
        graph[endWaypoint] = new Dictionary<GameObject, float>();

        if (!IsObstacleBetween(startWaypoint, endWaypoint))
        {
            float distance = Vector3.Distance(startWaypoint.transform.position, endWaypoint.transform.position);
            graph[startWaypoint][endWaypoint] = distance;
            graph[endWaypoint][startWaypoint] = distance;
            Debug.Log("Added edge from " + startWaypoint.name + " to " + endWaypoint.name + " with distance " + distance);
        }
    }

    bool IsObstacleBetween(GameObject waypoint1, GameObject waypoint2)
    {
        Vector3 direction = waypoint2.transform.position - waypoint1.transform.position;
        float distance = direction.magnitude;
        Ray ray = new Ray(waypoint1.transform.position, direction);
        return Physics.Raycast(ray, distance);
    }

    void FindPath(GameObject start, GameObject end)
    {
        Dictionary<GameObject, float> distances = new Dictionary<GameObject, float>();
        Dictionary<GameObject, GameObject> previous = new Dictionary<GameObject, GameObject>();
        List<GameObject> unvisited = new List<GameObject>(graph.Keys);

        foreach (var node in graph.Keys)
        {
            distances[node] = float.PositiveInfinity;
            previous[node] = null;
        }

        distances[start] = 0;

        while (unvisited.Count > 0)
        {
            GameObject currentNode = null;
            float smallestDistance = float.PositiveInfinity;

            foreach (var node in unvisited)
            {
                if (distances[node] < smallestDistance)
                {
                    smallestDistance = distances[node];
                    currentNode = node;
                }
            }

            if (currentNode == null)
                break;

            unvisited.Remove(currentNode);

            foreach (var neighbor in graph[currentNode].Keys)
            {
                float newDistance = distances[currentNode] + graph[currentNode][neighbor];

                if (newDistance < distances[neighbor])
                {
                    distances[neighbor] = newDistance;
                    previous[neighbor] = currentNode;
                }
            }
        }

        List<GameObject> path = new List<GameObject>();
        GameObject current = end;

        while (current != null)
        {
            path.Add(current);
            current = previous[current];
        }

        path.Reverse();
        Debug.Log("Path found: " + string.Join(" -> ", path.ConvertAll(node => node.name).ToArray()));
        VisualizePath(path);
    }

     void VisualizePath(List<GameObject> path)
    {
        Debug.DrawLine(startWaypoint.transform.position, endWaypoint.transform.position, Color.red, 25f);
        Debug.Log("Drawing line from " + startWaypoint.name + " to " + endWaypoint.name);

        /*for (int i = 0; i < path.Count - 1; i++)
        {
            Debug.DrawLine(startWaypoint.transform.position, endWaypoint.transform.position, Color.red, 25f);
            Debug.Log("Drawing line from " + startWaypoint.name + " to " + endWaypoint.name);
        }*/
    }
}

