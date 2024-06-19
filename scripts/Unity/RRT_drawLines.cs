using UnityEngine;
using System.Collections.Generic;

public class RRT_drawlines : MonoBehaviour
{
    //public Text statusText; // UI text to display status (optional)
    public GameObject startWaypoint; // Start waypoint assigned in the Inspector
    public GameObject endWaypoint; // End waypoint assigned in the Inspector
    public List<GameObject> obstacles; // List of all obstacle objects assigned in the Inspector
    private List<GameObject> newPoints = new List<GameObject>();
    private Dictionary<GameObject, Dictionary<GameObject, float>> graph = new Dictionary<GameObject, Dictionary<GameObject, float>>();

    void Start()
    {
        Debug.Log("Start: Checking assignments before initialization");
        Debug.Log("startWaypoint: " + (startWaypoint != null ? startWaypoint.name : "null"));
        Debug.Log("endWaypoint: " + (endWaypoint != null ? endWaypoint.name : "null"));

        AddMeshCollidersToObstacles();
        AddMeshCollidersToStartEndPoints();

        if (startWaypoint == null || endWaypoint == null)
        {
            Debug.LogError("Start or End waypoint is not assigned. Please assign waypoints in the Inspector.");
            return;
        }

        GenerateRRTPoints(startWaypoint, endWaypoint, 5, 4);
        Debug.Log("RRTPoints found.");

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

    void AddMeshCollidersToStartEndPoints()
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

    void AddMeshColliderToPoint(GameObject point)
    {
        MeshCollider meshCollider = point.GetComponent<MeshCollider>();
        if (meshCollider == null)
        {
            meshCollider = point.AddComponent<MeshCollider>();
            Debug.Log("Mesh Collider added to " + point.name);
        }
    }

    bool IsObstacleBetween(GameObject waypoint1, GameObject waypoint2)
    {
        Vector3 direction = waypoint2.transform.position - waypoint1.transform.position;
        float distance = direction.magnitude;
        Ray ray = new Ray(waypoint1.transform.position, direction);
        int layerMask = LayerMask.GetMask("Obstacles"); // Ensure obstacles are on this layer

        return Physics.Raycast(ray, distance, layerMask);
    }

    void GenerateRRTPoints(GameObject start, GameObject target, int numberOfPoints, float maxDistance)
    {
        GameObject currentPoint = start;

        for (int i = 0; i < numberOfPoints; i++)
        {
            Vector3 randomPoint = new Vector3(
                Random.Range(currentPoint.transform.position.x - maxDistance, currentPoint.transform.position.x + maxDistance),
                Random.Range(currentPoint.transform.position.y - maxDistance, currentPoint.transform.position.y + maxDistance),
                Random.Range(currentPoint.transform.position.z - maxDistance, currentPoint.transform.position.z + maxDistance)
            );

            GameObject newPoint = new GameObject("New Waypoint" + i);
            newPoint.transform.position = randomPoint;

            AddMeshColliderToPoint(newPoint);

            if (!IsObstacleBetween(currentPoint, newPoint))
            {
                float distanceToCurrent = Vector3.Distance(currentPoint.transform.position, newPoint.transform.position);

                if (!graph.ContainsKey(currentPoint))
                {
                    graph[currentPoint] = new Dictionary<GameObject, float>();
                }
                if (!graph.ContainsKey(newPoint))
                {
                    graph[newPoint] = new Dictionary<GameObject, float>();
                }

                graph[currentPoint][newPoint] = distanceToCurrent;
                graph[newPoint][currentPoint] = distanceToCurrent;

                currentPoint = newPoint;
                newPoints.Add(newPoint);
                Debug.Log("Added new waypoint " + newPoint.name + " closer to target.");
            }
            else
            {
                Destroy(newPoint);
                Debug.Log("New waypoint " + i + " is not valid due to obstacles.");
            }

            // Check if the current point is close enough to the target
            if (Vector3.Distance(currentPoint.transform.position, target.transform.position) < maxDistance)
            {
                float distanceToTarget = Vector3.Distance(currentPoint.transform.position, target.transform.position);
                graph[currentPoint][target] = distanceToTarget;
                graph[target][currentPoint] = distanceToTarget;
                Debug.Log("Connected to target waypoint.");
                break;
            }
        }
    }

    void InitializeGraph()
    {
        Debug.Log("Initializing graph...");
        graph[startWaypoint] = new Dictionary<GameObject, float>();
        graph[endWaypoint] = new Dictionary<GameObject, float>();

        foreach (var point in newPoints)
        {
            if (!graph.ContainsKey(point))
            {
                graph[point] = new Dictionary<GameObject, float>();
            }
        }

        // Connect startWaypoint to the first new point
        if (newPoints.Count > 0 && !IsObstacleBetween(startWaypoint, newPoints[0]))
        {
            float distance = Vector3.Distance(startWaypoint.transform.position, newPoints[0].transform.position);
            graph[startWaypoint][newPoints[0]] = distance;
            graph[newPoints[0]][startWaypoint] = distance;
            Debug.Log("Added edge from " + startWaypoint.name + " to " + newPoints[0].name + " with distance " + distance);
        }

        // Connect all new points to each other sequentially
        for (int i = 0; i < newPoints.Count - 1; i++)
        {
            if (!IsObstacleBetween(newPoints[i], newPoints[i + 1]))
            {
                float distance = Vector3.Distance(newPoints[i].transform.position, newPoints[i + 1].transform.position);
                graph[newPoints[i]][newPoints[i + 1]] = distance;
                graph[newPoints[i + 1]][newPoints[i]] = distance;
                Debug.Log("Added edge from " + newPoints[i].name + " to " + newPoints[i + 1].name + " with distance " + distance);
            }
        }

        // Connect the last new point to endWaypoint
        if (newPoints.Count > 0 && !IsObstacleBetween(newPoints[newPoints.Count - 1], endWaypoint))
        {
            float distance = Vector3.Distance(newPoints[newPoints.Count - 1].transform.position, endWaypoint.transform.position);
            graph[newPoints[newPoints.Count - 1]][endWaypoint] = distance;
            graph[endWaypoint][newPoints[newPoints.Count - 1]] = distance;
            Debug.Log("Added edge from " + newPoints[newPoints.Count - 1].name + " to " + endWaypoint.name + " with distance " + distance);
        }
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
        for (int i = 0; i < path.Count - 1; i++)
        {
            if (!IsObstacleBetween(path[i], path[i + 1]))
            {
                Debug.DrawLine(path[i].transform.position, path[i + 1].transform.position, Color.red, 25f);
                Debug.Log("Drawing line from " + path[i].name + " to " + path[i + 1].name);
            }
            else
            {
                Debug.Log("Obstacle detected between " + path[i].name + " and " + path[i + 1].name + ". Path not visualized.");
                return; // Stop visualization if an obstacle is detected
            }
        }
    }


}

